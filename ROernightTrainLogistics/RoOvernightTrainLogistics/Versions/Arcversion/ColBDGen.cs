using System;
using System.IO;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ILOG.Concert;
using ILOG.CPLEX;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.Distributions;
using NumSharp;
using NumSharp.Extensions;
using System.Diagnostics;

namespace RoOvernightTrainLogistics
{
    class DataStructure
    {
        #region
        public double[,] x_ij_val;
        public double righthandsideXi;
        public double[,] v_ij_val;
        public double[] sigma_j_val;
        public double[] w_i_val;
        public double[,] MXi_val;
        public DataStructure(Data data)
        {
            x_ij_val = new double[data.nodeSize, data.DCSize];
            v_ij_val = new double[data.nodeSize, data.DCSize];
            MXi_val = new double[data.nodeSize, data.DCSize];
            sigma_j_val = new double[data.DCSize];
            w_i_val = new double[data.nodeSize];
        }
        #endregion
    }
    /// <summary>
    /// embedding MILP for subproblem into B&C framework
    /// invoked at node which finds an integer feasible solution
    /// </summary>
    internal class SCNRLazyConsCallback : Cplex.LazyConstraintCallback
    {
        internal  SCNRMaster Scenariomaster;
        internal Data data;
        internal int[] y_solution;        
        internal SCNRSub modifiedSub;

        internal List<int[]> solutionpool;
        internal List<int> searchtimes;
        internal double bestsubcost;
        internal double bestmastercost;

        public SCNRLazyConsCallback(SCNRMaster Scenariomaster, Data data, int[] y_solution)
        {
            this.Scenariomaster = Scenariomaster;
            this.data = data;
            this.y_solution = y_solution;
            
            bestsubcost = -1;
            bestmastercost = -1;
            modifiedSub = new SCNRSub(y_solution, new int[data.linkSize], data);
            modifiedSub.GenScenarioSubproblem();
            solutionpool = new List<int[]>();
            searchtimes = new List<int>();

        }
        public override void Main()
        {
            #region
            int[] u_sol = new int[data.sectionSize];
            int[] u_sol_to_link_state = new int[data.linkSize];

            int numberofdisruptions = 0;
            for (int a = 0; a < data.sectionSize; a++)
            {
                double val = GetValue(Scenariomaster.u_a[a]);
                if (val > 0.5)
                {
                    u_sol[a] = 1;
                    numberofdisruptions++;
                }
            }
            for (int s = 0; s < data.sectionSize; s++)
            {
                //transfer section state to link state
                if (u_sol[s] == 1)
                {
                    Section section = data.sectionlist[s];
                    for (int a = 0; a < section.CorrespondLinks.Count; a++)
                    {
                        u_sol_to_link_state[section.CorrespondLinks[a].LinkSeqNo] = 1;
                    }
                }
            }
            //reset the objective of subproblem and the right hand side of corresponding constraints
            modifiedSub.ResetOBJaRHS_XI(y_solution, u_sol_to_link_state, data.LB_Vi);
            bool issolved = modifiedSub.model.Solve();
            modifiedSub.model.ExportModel("modifiedSub.lp");
            Cplex.Status solvingstatus = modifiedSub.model.GetStatus();

            if (issolved)
            {
                double[] x_r_values = new double[data.pathSize];
                bool isfeaisble = true;
                for (int r = 0; r < data.pathSize; r++)
                {
                    x_r_values[r] = modifiedSub.model.GetValue(modifiedSub.x_r[r]);
                    for (int a = 0; a < data.pathlist[r].Count; a++)
                    {
                        int linkno = data.pathlist[r][a];

                        if (u_sol_to_link_state[linkno] * x_r_values[r] >= data.epsilon)
                        {
                            Program.Solvinglog.WriteLine("Error during finding the scenario using lazyconstraint callback: " + x_r_values[r]);
                            isfeaisble = false;
                        }
                    }
                }
                if (isfeaisble)
                {
                    double benders_primal_cost = modifiedSub.model.GetObjValue();

                    if (benders_primal_cost > Scenariomaster.bestObjVal)
                    {
                        Scenariomaster.bestObjVal = benders_primal_cost;
                        Scenariomaster.bestFeasibleSolution = u_sol;
                    }
                    else
                    {
                        bool sametobestsolution = true;
                        for (int a = 0; a < data.sectionSize; a++)
                        {
                            if (Scenariomaster.bestFeasibleSolution[a] != u_sol[a])
                            {
                                sametobestsolution = false; break;
                            }
                        }

                        if (data.nogoodcuts && benders_primal_cost < Scenariomaster.bestObjVal)
                        {
                            ILinearNumExpr combinatorialCut = Scenariomaster.model.LinearNumExpr();

                            double cardi = 0;
                            for (int a = 0; a < data.sectionSize; a++)
                            {
                                if (u_sol[a] == 1)
                                {
                                    combinatorialCut.AddTerm(-1, Scenariomaster.u_a[a]);
                                    cardi++;
                                }
                                else
                                {
                                    combinatorialCut.AddTerm(1, Scenariomaster.u_a[a]);
                                }
                            }
                            Add(Scenariomaster.model.Ge(combinatorialCut, 1 - cardi));

                        }
                    }

                    ILinearNumExpr extrem_point_exp = Scenariomaster.model.LinearNumExpr();

                    for (int r = 0; r < data.pathSize; r++)
                    {
                        for (int a = 0; a < data.pathlist[r].Count; a++)
                        {
                            int linkno = data.pathlist[r][a];
                            int sectionno = Program.g_link_list[linkno].CoreTSectionNo;
                            extrem_point_exp.AddTerm(-data.LB_Vi[data.e_r[r]] * x_r_values[r], Scenariomaster.u_a[sectionno]);
                        }

                    }
                    //add terms to objective                   
                    extrem_point_exp.AddTerm(1, Scenariomaster.eta);
                    //add benders cut
                    Add(Scenariomaster.model.Le(extrem_point_exp, benders_primal_cost));

                }
                else
                {
                    ILinearNumExpr combinatorialCut = Scenariomaster.model.LinearNumExpr();

                    double cardi = 0;
                    for (int a = 0; a < data.sectionSize; a++)
                    {
                        if (u_sol[a] == 1)
                        {
                            combinatorialCut.AddTerm(-1, Scenariomaster.u_a[a]);
                            cardi++;
                        }
                        else
                        {
                            combinatorialCut.AddTerm(1, Scenariomaster.u_a[a]);
                        }
                    }
                    Add(Scenariomaster.model.Ge(combinatorialCut, 1 - cardi));
                }
            }
            
            #endregion
        }
    }
    class SCNRMaster
    {
        #region
        public Cplex model;

        public INumVar eta;
        public INumVar[] u_a;
        Data data = new Data();
        public int number_of_var;
        public int number_of_con;

        public List<IRange> cutsStore;
        public List<IRange> optimalityCuts;
        public List<DataStructure> XiCoeff_list;
        public List<DataStructure> Temp_XiCoeff_list;
        //used for checking the correctness of B&C
        public int[] bestFeasibleSolution;
        public double bestObjVal = -float.MaxValue;
        public List<double> objvalpool = new List<double>();
        public List<int[]> solutionpool = new List<int[]>();

        internal int[] y_sol;
        public SCNRMaster(int[] y_sol)
        {
            this.y_sol = y_sol;
        }
        /// <summary>
        /// generate MILP model for scenario master problem
        /// </summary>
        /// <param name="y_sol"></param>
        public void GenScenarioMaster()
        {
            #region
            model = new Cplex();

            u_a = new INumVar[data.sectionSize];

            optimalityCuts = new List<IRange>();
            cutsStore = new List<IRange>();

            eta = model.NumVar(-float.MaxValue, float.MaxValue, NumVarType.Float);

            for (int a = 0; a < data.sectionSize; a++)
            {
                // set the upperbound of xi_j as the value of y_j
                u_a[a] = model.BoolVar($"u_{a}");
                number_of_var++;
            }
            model.AddMaximize(eta);

            int cons_1 = 1;
            if (cons_1 == 1)
            {
                ILinearNumExpr constraint = model.LinearNumExpr();
                ILinearNumExpr constraint1 = model.LinearNumExpr();
                int[] sectionstate = new int[data.sectionSize];

                for (int r = 0; r < data.pathSize; r++)
                {
                    if(y_sol[data.s_r[r]] == 1)
                    {
                        for (int a = 0; a < data.pathlist[r].Count; a++)
                        {
                            int linkno = data.pathlist[r][a];
                            int sectionno = Program.g_link_list[linkno].CoreTSectionNo;
                            sectionstate[sectionno] = 1;
                        }
                    }
                }

                for (int a = 0; a < data.sectionSize; a++)
                {
                    if(sectionstate[a] == 1)
                    {
                        constraint.AddTerm(1, u_a[a]);
                    }
                    else
                    {
                        constraint1.AddTerm(1, u_a[a]);
                    }
                }

                model.AddEq(constraint, data.max_dstroyed_DCs, "Disruption budget");
                model.AddEq(constraint1, 0, "Disruption budget");
                number_of_con++;
            }

            #endregion
        }
        /// <summary>
        /// generate linear model for scenario master problem
        /// \xi is continuous here!!
        /// to stabilize a new linear \xi 
        /// </summary>
        /// <param name="y_sol"></param>
        public void GenInternalPoint()
        {
            #region
            model = new Cplex();

            u_a = new INumVar[data.linkSize];

            optimalityCuts = new List<IRange>();
            cutsStore = new List<IRange>();
            XiCoeff_list = new List<DataStructure>();
            for (int a = 0; a < data.linkSize; a++)
            {
                // set the upperbound of xi_j as the value of y_j
                u_a[a] = model.BoolVar($"u_{a}");
                number_of_var++;
            }

            ILinearNumExpr objexpr = model.LinearNumExpr();
            for (int a = 0; a < data.linkSize; a++)
            {
                objexpr.AddTerm(1, u_a[a]);

                number_of_var++;
            }
            model.AddMaximize(objexpr);

            int cons_1 = 1;
            if (cons_1 == 1)
            {
                ILinearNumExpr constraint = model.LinearNumExpr();

                for (int a = 0; a < data.linkSize; a++)
                {
                    constraint.AddTerm(1, u_a[a]);
                }

                model.AddEq(constraint, data.max_dstroyed_DCs, "Disruption budget");
                number_of_con++;
            }
            #endregion
        }
        /// <summary>
        /// check the feasibility of a solution to judge if a no-good cut is added
        /// </summary>
        /// <param name="u_sol"></param>
        public bool FeasibilityCheck(int[]y_sol, int[] u_sol_to_link_state)
        {
            #region
            bool feasi = true;
            int[] nodedemandstatus = new int[data.nodeSize];
            Dictionary<int, List<int>> demandnodeno_to_servicepathlist_dic = new Dictionary<int, List<int>>();
            for (int r = 0; r < data.pathSize; r++)
            {
                if (y_sol[data.s_r[r]] == 1)
                {
                    nodedemandstatus[data.e_r[r]] += 1;
                    if (demandnodeno_to_servicepathlist_dic.ContainsKey(data.e_r[r]))
                    {
                        demandnodeno_to_servicepathlist_dic[data.e_r[r]].Add(r);
                    }
                    else
                    {
                        List<int> temp_list = new List<int>();
                        temp_list.Add(r);
                        demandnodeno_to_servicepathlist_dic.Add(data.e_r[r], temp_list);
                    }
                }
            }
            foreach(var item in demandnodeno_to_servicepathlist_dic)
            {
                
                //means only one path is used to serve the demandnode
                //then if one of the arc in the path is destroyed, then the solution is infeasible
                int numberofinfeasiroute = 0;
                for (int i = 0; i < item.Value.Count; i++)
                {
                    int r = item.Value[i];
                    
                    for (int a = 0; a < data.pathlist[r].Count; a++)
                    {
                        int linkno = data.pathlist[r][a];
                        if (u_sol_to_link_state[linkno] == 1)
                        {
                            numberofinfeasiroute++;                            
                        }
                    }
                }
                //if the number of destroyted path equals the size of available
                if (numberofinfeasiroute >= item.Value.Count)
                {
                    feasi = false; break;
                }                    
            }
            
            return feasi;
            #endregion
        }
        #endregion
    }
    class SCNRSub
    {
        #region
        public Cplex model;
        public Data data;

        internal INumVar[] x_r;
        internal INumVar[] z_j;
        
        public IRange[] relatedtodual_V_value;
        public IRange[] relatedtodual_W_value;
        public IRange[] relatedtodual_Delta_value;
        public IRange[] relatedtodual_Gamma_value;

        public int number_of_var;
        public int number_of_con;

        internal int[] u_sol_to_link_state;
        internal int[] y_sol;
        
        /// <param name="falocsolution"></param>        
        /// <param name="scenario_xi_j"></param>
        public SCNRSub(int[] y_sol, int[] u_sol_to_link_state, Data data)
        {
            this.y_sol = y_sol;            
            this.u_sol_to_link_state = u_sol_to_link_state;
            this.data = data;
        }
        /// <summary>
        /// ccg subproblem 
        /// used for the valid cut for ccg master problem
        /// </summary>        
        public void GenScenarioSubproblem()
        {
            #region ccg subproblem model           
            //z_j_sub_vars_arr = new INumVar[data.DCSize];

            model = new Cplex();            
                        
            //对于小规模问题可加速,大规模系数问题效率不高,实测可在25*15*4算例中加速
            //model.SetParam(Cplex.Param.Barrier.Crossover, 2);
            number_of_var = 0; number_of_con = 0;
            
            x_r = new INumVar[data.pathSize];
            z_j = new INumVar[data.DCSize];
            number_of_var = 0; number_of_con = 0;

            relatedtodual_V_value = new IRange[data.nodeSize];
            relatedtodual_W_value = new IRange[data.DCSize];
            relatedtodual_Delta_value = new IRange[data.pathSize];
            relatedtodual_Gamma_value = new IRange[data.DCSize];

            for (int r = 0; r < data.pathSize; r++)
            {
                x_r[r] = model.NumVar(0, float.MaxValue, NumVarType.Float, $"x_{r}");
                number_of_var++;
            }

            for (int j = 0; j < data.DCSize; j++)
            {
                z_j[j] = model.NumVar(0, float.MaxValue, NumVarType.Float, $"z_{j}");
                number_of_var++;                
            }
            
            ILinearNumExpr singleterm = model.LinearNumExpr();
            for (int j = 0; j < data.DCSize; j++)
            {                
                singleterm.AddTerm(data.c_j[j], z_j[j]);
                
            }
            for (int r = 0; r < data.pathSize; r++)
            {
                if(data.pathlist[r].Count != 0)
                {
                    for (int a = 0; a < data.pathlist[r].Count; a++)
                    {
                        int linkno = data.pathlist[r][a];
                        singleterm.AddTerm((data.d_r[r] + data.LB_Vi[data.e_r[r]] * u_sol_to_link_state[linkno]), x_r[r]);
                    }
                }
                else
                {
                    singleterm.AddTerm(data.d_r[r], x_r[r]);
                }             
            }
            //add terms to objective
            model.AddMinimize(singleterm);

            //constraints      
            
            //constraints            
            int cons_1 = 1;
            if (cons_1 == 1)
            {
                for (int r = 0; r < data.pathSize; r++)
                {                    
                    IRange constraint = model.AddLe(x_r[r], y_sol[data.s_r[r]], $"Transport_capacity_{r}");

                    relatedtodual_Delta_value[r] = constraint;                    
                }
            }

            int cons_2 = 1;
            if (cons_2 == 1)
            {
                for (int i = 0; i < data.nodeSize; i++)
                {
                    ILinearNumExpr contraint = model.LinearNumExpr();
                    for (int r = 0; r < data.pathSize; r++)
                    {                        
                        if(data.e_r[r] == i)
                        {
                            contraint.AddTerm(1, x_r[r]);
                        }
                    }
                    
                    IRange constraint = model.AddGe(contraint, 1, $"2st Fulfill_demand_{i}");
                    relatedtodual_V_value[i] = constraint;
                    number_of_con++;
                }
            }            
            int cons_3 = 1;
            if (cons_3 == 1)
            {
                for (int j = 0; j < data.DCSize; j++)
                {
                    ILinearNumExpr contraint = model.LinearNumExpr();
                    for (int r = 0; r < data.pathSize; r++)
                    {                        
                        if (data.s_r[r] > data.DCSize)
                        {
                            Console.WriteLine("Find an error here!!!!!");
                        }

                        if (data.s_r[r] == j)
                        {
                            contraint.AddTerm(data.h_i[data.e_r[r]], x_r[r]);
                        }
                    }
                    contraint.AddTerm(-1, z_j[j]);
                    IRange constraint = model.AddLe(contraint,0, $"purchasinggoods_{j}");
                    relatedtodual_W_value[j] = constraint;
                    number_of_con++;
                }
            }
            int cons_4 = 1;
            if (cons_4 == 1)
            {
                int[] openstate = new int[data.DCSize];

                for (int a = 0; a < data.linkSize; a++)
                {
                    if(u_sol_to_link_state[a] == 1)
                    {
                        int fromnode = Program.g_link_list[a].InternalFromNodeSeqNo;
                        int tonode = Program.g_link_list[a].InternalToNodeSeqNo;

                        int dcnodefrom = Program.g_node_list[fromnode].CorDcenterSeqNo;
                        int dcnodeto = Program.g_node_list[tonode].CorDcenterSeqNo;

                        openstate[dcnodefrom] = 1;
                        openstate[dcnodeto] = 1;
                    }
                }
                
                for (int j = 0; j < data.DCSize; j++)
                {
                    IRange constraint = model.AddLe(z_j[j], y_sol[j] * (1 - openstate[j]) * data.big_M, $"purchasinggoods_z_{j}");
                    relatedtodual_Gamma_value[j] = constraint;
                    number_of_con++;
                }
            }

            #endregion
        }
        /// <summary>
        /// Update model at each branch or iteration without newly formulating
        /// </summary>
        /// <param name="upd_y_solution"></param-difference>        
        /// <param name="upd_u_sol"></param-difference>
        public void ResetOBJaRHS_XI(int[] upd_y_solution, int[] upd_u_sol_to_link_state, double[] LB_Vi)
        {
            #region

            IObjective mdsubobj = model.GetObjective();

            ILinearNumExpr singleterm = model.LinearNumExpr();
            for (int j = 0; j < data.DCSize; j++)
            {
                singleterm.AddTerm(data.c_j[j], z_j[j]);

            }
            for (int r = 0; r < data.pathSize; r++)
            {
                if (data.pathlist[r].Count != 0)
                {
                    for (int a = 0; a < data.pathlist[r].Count; a++)
                    {
                        int linkno = data.pathlist[r][a];
                        singleterm.AddTerm((data.d_r[r] + LB_Vi[data.e_r[r]] * upd_u_sol_to_link_state[linkno]), x_r[r]);
                    }
                }
                else
                {
                    singleterm.AddTerm(data.d_r[r], x_r[r]);
                }
            }            
            //clear the current objective
            mdsubobj.ClearExpr();
            //reset the objective
            mdsubobj.Expr = singleterm;

            //Reset right hand side of constraints related to y z xi solutions
            for (int r = 0; r < data.pathSize; r++)
            {                
                if (relatedtodual_Delta_value[r] != null)
                {
                    relatedtodual_Delta_value[r].UB = upd_y_solution[data.s_r[r]];
                }           
            }

            int[] openstate = new int[data.DCSize];

            for (int a = 0; a < data.linkSize; a++)
            {
                if (upd_u_sol_to_link_state[a] == 1)
                {
                    int fromnode = Program.g_link_list[a].InternalFromNodeSeqNo;
                    int tonode = Program.g_link_list[a].InternalToNodeSeqNo;

                    int dcnodefrom = Program.g_node_list[fromnode].CorDcenterSeqNo;
                    int dcnodeto = Program.g_node_list[tonode].CorDcenterSeqNo;

                    openstate[dcnodefrom] = 1;
                    openstate[dcnodeto] = 1;
                }
            }

            for (int j = 0; j < data.DCSize; j++)
            {
                relatedtodual_Gamma_value[j].UB = y_sol[j] * (1 - openstate[j]) * data.big_M;
                number_of_con++;
            }
            #endregion
        }
        #endregion
    }
    class CCGMaster
    {
        public Cplex model;

        public INumVar omega;
        public INumVar[] y_j;
        
        public List<INumVar[]> x_l_r;
        public List<INumVar[]> z_l_j;
        public List<int[]> u_l_a_to_link_state;

        public int number_of_var;
        public int number_of_con;

        List<IRange> generateCuts;
        public List<DataStructure> YCoeff_list;
        public List<DataStructure> Temp_YCoeff_list;

        Data data = new Data();
        public void InitializeCCGMaster()
        {
            #region model
            model = new Cplex();

            number_of_var = 0; number_of_con = 0;

            y_j = new INumVar[data.DCSize];
            
            x_l_r = new List<INumVar[]>();
            z_l_j = new List<INumVar[]>();
            u_l_a_to_link_state = new List<int[]>();

            for (int j = 0; j < data.DCSize; j++)
            {
                y_j[j] = model.BoolVar($"y_{j}");
                number_of_var++;
            }

            omega = model.NumVar(float.MinValue, float.MaxValue, NumVarType.Float, "omega");

            ILinearNumExpr numExpr = model.LinearNumExpr();
            for (int j = 0; j < data.DCSize; j++)
            {
                numExpr.AddTerm(data.f_j[j], y_j[j]);                
            }
            numExpr.AddTerm(1, omega);            
            model.AddMinimize(numExpr);
            
            int cons_1 = 0;
            if (cons_1 == 1)
            {
                ILinearNumExpr totalOpenDCs = model.LinearNumExpr();
                for (int j = 0; j < data.DCSize; j++)
                {
                    totalOpenDCs.AddTerm(1, y_j[j]);
                }
                model.AddLe(totalOpenDCs, data.max_open_DCs, $"open size");
            }
            int cons_2 = 1;
            if (cons_2 == 1)
            {
                ILinearNumExpr totalOpenDCs = model.LinearNumExpr();
                for (int j = 0; j < data.DCSize; j++)
                {
                    totalOpenDCs.AddTerm(1, y_j[j]);
                }
                model.AddGe(totalOpenDCs, data.Min_open_DCs, $"leastDCs");
            }
            #endregion
        }
        /// <summary>
        /// ccg master problem
        /// </summary>
        /// <param name="cutting_indicator"></param>
        public void GenCCGMasterproblem(int cutting_indicator)
        {
            #region model
            //add the newly finded scenario into CCG model
            int[] worstscenario = u_l_a_to_link_state.Last();

            number_of_var = 0; number_of_con = 0;

            INumVar[] temp_x_r = new INumVar[data.pathSize];
            INumVar[] temp_z_j = new INumVar[data.DCSize];
            number_of_var = 0; number_of_con = 0;

            for (int r = 0; r < data.pathSize; r++)
            {
                temp_x_r[r] = model.NumVar(0, float.MaxValue, NumVarType.Float, $"x_{x_l_r.Count}_{r}");
                number_of_var++;
            }
            x_l_r.Add(temp_x_r);

            for (int j = 0; j < data.DCSize; j++)
            {
                temp_z_j[j] = model.NumVar(0, float.MaxValue, NumVarType.Float, $"z_{x_l_r.Count}_{j}");
                number_of_var++;
            }
            z_l_j.Add(temp_z_j);

            if(cutting_indicator == 1)
            {
                ILinearNumExpr singleterm = model.LinearNumExpr();
                for (int j = 0; j < data.DCSize; j++)
                {
                    singleterm.AddTerm(data.c_j[j], temp_z_j[j]);

                }
                for (int r = 0; r < data.pathSize; r++)
                {
                    singleterm.AddTerm(data.d_r[r], temp_x_r[r]);
                }
                //add terms to objective
                model.AddGe(omega, singleterm, "optimality cut");
            }

            //constraints      
           
            //constraints            
            int cons_1 = 1;
            if (cons_1 == 1)
            {
                for (int r = 0; r < data.pathSize; r++)
                {
                    if(data.pathlist[r].Count != 0)
                    {
                        for (int a = 0; a < data.pathlist[r].Count; a++)
                        {
                            int linkno = data.pathlist[r][a];

                            model.AddLe(temp_x_r[r], model.Prod(y_j[data.s_r[r]], (1 - worstscenario[linkno])), $"Transport_capacity_{r}_{linkno}");

                        }
                    }
                    else
                    {
                        model.AddLe(temp_x_r[r], y_j[data.s_r[r]], $"Transport_capacity_{r}_null");
                    }
                }
            }

            int cons_2 = 1;
            if (cons_2 == 1)
            {
                for (int i = 0; i < data.nodeSize; i++)
                {
                    ILinearNumExpr contraint = model.LinearNumExpr();
                    for (int r = 0; r < data.pathSize; r++)
                    {
                        if (data.e_r[r] == i)
                        {
                            contraint.AddTerm(1, temp_x_r[r]);
                        }
                    }

                    IRange constraint = model.AddGe(contraint, 1, $"2st Fulfill_demand_{i}");
                    
                    number_of_con++;
                }
            }
            int cons_3 = 1;
            if (cons_3 == 1)
            {
                for (int j = 0; j < data.DCSize; j++)
                {
                    ILinearNumExpr contraint = model.LinearNumExpr();
                    for (int r = 0; r < data.pathSize; r++)
                    {                        
                        if (data.s_r[r] == j)
                        {
                            contraint.AddTerm(data.h_i[data.e_r[r]], temp_x_r[r]);
                        }
                    }
                    contraint.AddTerm(-1, temp_z_j[j]);
                    IRange constraint = model.AddLe(contraint, 0, $"purchasinggoods_{j}");
                    
                    number_of_con++;
                }
            }
            int cons_4 = 1;
            if (cons_4 == 1)
            {
                int[] openstate = new int[data.DCSize];

                for (int a = 0; a < data.linkSize; a++)
                {
                    if (worstscenario[a] == 1)
                    {
                        int fromnode = Program.g_link_list[a].InternalFromNodeSeqNo;
                        int tonode = Program.g_link_list[a].InternalToNodeSeqNo;

                        int dcnodefrom = Program.g_node_list[fromnode].CorDcenterSeqNo;
                        int dcnodeto = Program.g_node_list[tonode].CorDcenterSeqNo;

                        openstate[dcnodefrom] = 1;
                        openstate[dcnodeto] = 1;
                    }
                }

                for (int j = 0; j < data.DCSize; j++)
                {
                    model.AddLe(temp_z_j[j], model.Prod(y_j[j], (1 - openstate[j]) * data.big_M), $"purchasinggoods_z_{j}");
                    number_of_con++;
                }
            }
            #endregion
        }
        
    }
    class CCGSub
    {
        #region
        public Cplex model;
        public Data data = new Data();

        internal INumVar[] x_r;
        internal INumVar[] z_j;
        internal INumVar[] p_i;

        public IRange[] relatedtodual_V_value;
        public IRange[] relatedtodual_W_value;
        public IRange[] relatedtodual_Delta_value;

        public int number_of_var;
        public int number_of_con;
        
        /// <summary>
        /// ccg subproblem 
        /// used for the valid cut for ccg master problem
        /// </summary>        
        public void GenCCGSubproblem(int[] y_sol, int[] u_sol_to_link_state)
        {
            #region ccg subproblem model           
            
            model = new Cplex();
            //对于小规模问题可加速,大规模系数问题效率不高,实测可在25*15*4算例中加速
            //model.SetParam(Cplex.Param.Barrier.Crossover, 2);
           
            number_of_var = 0; number_of_con = 0;

            x_r = new INumVar[data.pathSize];
            z_j = new INumVar[data.DCSize];
            p_i = new INumVar[data.nodeSize];
            number_of_var = 0; number_of_con = 0;

            relatedtodual_V_value = new IRange[data.nodeSize];
            relatedtodual_W_value = new IRange[data.DCSize];
            relatedtodual_Delta_value = new IRange[data.pathSize];

            for (int r = 0; r < data.pathSize; r++)
            {
                x_r[r] = model.NumVar(0, float.MaxValue, NumVarType.Float, $"x_{r}");
                number_of_var++;
            }

            for (int j = 0; j < data.DCSize; j++)
            {
                z_j[j] = model.NumVar(0, float.MaxValue, NumVarType.Float, $"z_{j}");
                number_of_var++;
            }
            for (int i = 0; i < data.nodeSize; i++)
            {
                p_i[i] = model.NumVar(0, 1, NumVarType.Float, $"p_{i}");
                number_of_var++;
            }
            ILinearNumExpr singleterm = model.LinearNumExpr();
            for (int j = 0; j < data.DCSize; j++)
            {
                singleterm.AddTerm(data.c_j[j], z_j[j]);

            }
            for (int r = 0; r < data.pathSize; r++)
            {
                singleterm.AddTerm(data.d_r[r], x_r[r]);
            }
            for (int i = 0; i < data.nodeSize; i++)
            {
                singleterm.AddTerm((data.maxopencost + data.maxtransportcost + data.maxpurchasecost) * 5 * data.h_i[i], p_i[i]);
            }
            //add terms to objective
            model.AddMinimize(singleterm);

            //constraints         
            
            //constraints            
            int cons_1 = 1;
            if (cons_1 == 1)
            {
                for (int r = 0; r < data.pathSize; r++)
                {                   
                    if (data.pathlist[r].Count != 0)
                    {
                        for (int a = 0; a < data.pathlist[r].Count; a++)
                        {
                            int linkno = data.pathlist[r][a];

                            IRange constraint = model.AddLe(x_r[r], y_sol[data.s_r[r]] * (1 - u_sol_to_link_state[linkno]), $"Transport_capacity_{r}_{linkno}");

                            relatedtodual_Delta_value[r] = constraint;
                        }
                    }
                    else
                    {
                        model.AddLe(x_r[r], y_sol[data.s_r[r]], $"Transport_capacity_{r}_null");
                    }
                }
            }

            int cons_2 = 1;
            if (cons_2 == 1)
            {
                for (int i = 0; i < data.nodeSize; i++)
                {
                    ILinearNumExpr contraint = model.LinearNumExpr();
                    for (int r = 0; r < data.pathSize; r++)
                    {
                        if (data.e_r[r] == i)
                        {
                            contraint.AddTerm(1, x_r[r]);
                        }
                    }
                    contraint.AddTerm(1, p_i[i]);
                    IRange constraint = model.AddGe(contraint, 1, $"2st Fulfill_demand_{i}");
                    relatedtodual_V_value[i] = constraint;
                    number_of_con++;
                }
            }
            int cons_3 = 1;
            if (cons_3 == 1)
            {
                for (int j = 0; j < data.DCSize; j++)
                {
                    ILinearNumExpr contraint = model.LinearNumExpr();
                    for (int r = 0; r < data.pathSize; r++)
                    {
                        if (data.s_r[r] > data.DCSize)
                        {
                            Console.WriteLine("Find an error here!!!!!");
                        }

                        if (data.s_r[r] == j)
                        {
                            contraint.AddTerm(data.h_i[data.e_r[r]], x_r[r]);
                        }
                    }
                    contraint.AddTerm(-1, z_j[j]);
                    IRange constraint = model.AddLe(contraint, 0, $"purchasinggoods_{j}");
                    relatedtodual_W_value[j] = constraint;
                    number_of_con++;
                }
            }
            int cons_4 = 1;
            if (cons_4 == 1)
            {
                for (int r = 0; r < data.pathSize; r++)
                {                    
                    if (data.pathlist[r].Count != 0)
                    {
                        for (int a = 0; a < data.pathlist[r].Count; a++)
                        {
                            int linkno = data.pathlist[r][a];

                            model.AddLe(z_j[data.e_r[r]], y_sol[data.s_r[r]] * (1 - u_sol_to_link_state[linkno])*data.big_M, $"Transport_capacity_{r}_{linkno}");                            
                        }
                    }
                    else
                    {
                        model.AddLe(x_r[r], y_sol[data.s_r[r]], $"Transport_capacity_{r}_null");
                    }
                }
            }
            #endregion
        }
        /// <summary>
        /// Update model at each branch or iteration without newly formulating
        /// </summary>
        /// <param name="upd_y_solution"></param-difference>        
        /// <param name="upd_u_sol"></param-difference>
        public void ResetOBJaRHS_XI(int[] upd_y_solution, int[] upd_u_sol_to_link_state, double[] bigM)
        {
            #region

            IObjective mdsubobj = model.GetObjective();

            ILinearNumExpr singleterm = model.LinearNumExpr();
            for (int j = 0; j < data.DCSize; j++)
            {
                singleterm.AddTerm(data.c_j[j], z_j[j]);

            }

            for (int r = 0; r < data.pathSize; r++)
            {
                if (data.pathlist[r].Count != 0)
                {
                    for (int a = 0; a < data.pathlist[r].Count; a++)
                    {
                        int linkno = data.pathlist[r][a];
                        singleterm.AddTerm((data.d_r[r] + bigM[data.e_r[r]] * upd_u_sol_to_link_state[linkno]), x_r[r]);
                    }
                }
                else
                {
                    singleterm.AddTerm(data.d_r[r], x_r[r]);
                }                
            }

            //clear the current objective
            mdsubobj.ClearExpr();
            //reset the objective
            mdsubobj.Expr = singleterm;

            //Reset right hand side of constraints related to y z xi solutions
            for (int r = 0; r < data.pathSize; r++)
            {
                if (relatedtodual_Delta_value[r] != null)
                {
                    relatedtodual_Delta_value[r].UB = upd_y_solution[data.s_r[r]];
                }
            }
            #endregion
        }
        #endregion
    }
    class DualCCGSub
    {
        public Cplex model;
        public Data data = new Data();

        public INumVar[] v_i;        
        public INumVar[] w_j;
        public INumVar[,] delta_ra;        

        public int number_of_var;
        public int number_of_con;
        
        public List<DataStructure> XiCoeff_list;
        public List<DataStructure> Temp_XiCoeff_list;
        List<IRange> generateCuts;
        int[] y_solution;
        /// <summary>
        /// dual problem of the ccg subproblem
        /// used for generating the worst scenario \xi_j
        /// note that this dual subproblem cannot provide valid cuts
        /// </summary>
        public DualCCGSub(int[] y_solution)
        {
            this.y_solution = y_solution;
        }
        public void GenDualofCCGsub()
        {
            #region
            model = new Cplex();
            number_of_var = 0; number_of_con = 0;

            model.SetParam(Cplex.Param.RootAlgorithm, Cplex.Algorithm.Primal);
            
            v_i = new INumVar[data.nodeSize];
            w_j = new INumVar[data.DCSize];            
            delta_ra = new INumVar[data.pathSize, data.linkSize];
            
            for (int r = 0; r < data.pathSize; r++)
            {
                for (int a = 0; a < data.linkSize; a++)
                {
                    delta_ra[r, a] = model.NumVar(0, float.MaxValue, NumVarType.Float, $"delta_{r}_{a}");
                    number_of_var++;
                }
            }
            for (int i = 0; i < data.nodeSize; i++)
            {
                v_i[i] = model.NumVar(0, float.MaxValue, NumVarType.Float, $"v_{i}");
                number_of_var++;                
            }
            for (int j = 0; j < data.DCSize; j++)
            {
                w_j[j] = model.NumVar(float.MinValue, 0, NumVarType.Float, $"v_{j}");
                number_of_var++;
            }
            
            //first part            
            ILinearNumExpr subobj = model.LinearNumExpr();
            for (int r = 0; r < data.pathSize; r++)
            {
                for (int a = 0; a < data.pathlist[r].Count; a++)
                {
                    int linkno = data.pathlist[r][a];
                    subobj.AddTerm(0, delta_ra[r, linkno]);                    
                }                
            }

            for (int i = 0; i < data.nodeSize; i++)
            {
                subobj.AddTerm(1, v_i[i]);
            }

            model.AddMaximize(subobj);

            int cons_0 = 1;
            if (cons_0 == 1)
            {
                ILinearNumExpr consexpr = model.LinearNumExpr();
                for (int r = 0; r < data.pathSize; r++)
                {
                    consexpr.AddTerm(data.h_i[data.s_r[r]], w_j[data.s_r[r]]);
                    consexpr.AddTerm(1, v_i[data.e_r[r]]);

                    for (int a = 0; a < data.pathlist[r].Count; a++)
                    {
                        int linkno = data.pathlist[r][a];
                        consexpr.AddTerm(-1, delta_ra[r, linkno]);
                    }
                    model.AddLe(consexpr, data.d_r[r], $"X_dual[{r}]");
                }
            }

            int cons_1 = 1;
            if (cons_1 == 1)
            {
                for (int j = 0; j < data.DCSize; j++)
                {                    
                    IRange constraint = model.AddGe(w_j[j], -data.c_j[j], $"Z_Dual_{j}");
                    number_of_con++;
                }
            }      
            #endregion
        }        
    }    
    class CACGBD
    {
        static Data data = new Data();
        Solution solution = new Solution();

        public List<int[]> xI_sols_pool;
        public double scenarioCost;

        /// <summary>
        /// transfer the current to a new model for the consideration 
        /// of branch and check scheme when facing a new schenario
        /// </summary>
        public SCNRMaster BranchandCutForWorstCsenario(int[] y_sol)
        {
            #region
            
            SCNRMaster ssmaster = new SCNRMaster(y_sol);
            xI_sols_pool = new List<int[]>();//save suboptimal scenarios
            int[] usol_to_link_state = new int[data.linkSize];
            
            for (int r = 0; r < data.pathSize; r++)
            {
                for (int a = 0; a < data.pathlist[r].Count; a++)
                {
                    int linkno = data.pathlist[r][a];
                    usol_to_link_state[linkno] = 1;
                }
            }
            CCGSub ccgsub = new CCGSub();
            ccgsub.GenCCGSubproblem(y_sol, usol_to_link_state);
            ccgsub.model.SetParam(Cplex.Param.Preprocessing.Presolve, false);
            //ccgsub.model.SetParam(Cplex.Param.RootAlgorithm, Cplex.Algorithm.Dual);
            ccgsub.model.Solve();

            ccgsub.model.ExportModel("CCGsubproblem.lp");

            double[] xvalue = new double[data.pathSize];
            for (int r = 0; r < data.pathSize; r++)
            {
                xvalue[r] = ccgsub.model.GetValue(ccgsub.x_r[r]);
            }
            for (int j = 0; j < data.DCSize; j++)
            {
                Console.WriteLine(ccgsub.model.GetValue(ccgsub.z_j[j]));
            }
            Cplex.Status solvingstatus = ccgsub.model.GetStatus();
            
            for (int i = 0; i < data.nodeSize; i++)
            {
                data.LB_Vi[i] = ccgsub.model.GetDual(ccgsub.relatedtodual_V_value[i]);
            }

            //running the final MILP model
            ssmaster.GenScenarioMaster();//generate master problem

            if (data.cplex_para_set)
            {
                // 设置全强分支策略
                ssmaster.model.SetParam(Cplex.Param.MIP.Strategy.VariableSelect, 3); // CPX_PARAM_VARSEL=3

                // 设置纯最佳边界搜索策略
                ssmaster.model.SetParam(Cplex.Param.MIP.Strategy.NodeSelect, 1); // CPX_PARAM_NODESEL=1
                ssmaster.model.SetParam(Cplex.Param.MIP.Interval, 1); // CPX_PARAM_BBINTERVAL=1
                ssmaster.model.SetParam(Cplex.Param.MIP.Tolerances.UpperCutoff, 0.0); // CPX_PARAM_BTTOL=0.0

                // 设置RINS启发式在每个节点执行
                ssmaster.model.SetParam(Cplex.Param.MIP.Limits.RepairTries, 1); // CPX_PARAM_RINSHEUR=1
            }

            // Turn on traditional search for use with control callbacks
            ssmaster.model.SetParam(Cplex.Param.MIP.Strategy.Search, Cplex.MIPSearch.Traditional);
            ssmaster.model.SetParam(Cplex.Param.Preprocessing.Presolve, false);
            ssmaster.model.SetParam(Cplex.Param.Threads, 1);

            //necessary, otherwise, GetObjValue differs from GetBestObjValue
            //ssmaster.model.SetParam(Cplex.Param.MIP.Tolerances.MIPGap, 1e-9);

            //get multicuts from solution pool
            ssmaster.model.SetParam(Cplex.Param.MIP.Pool.Capacity, data.maxMultiCuts); // 解池容量
            ssmaster.model.SetParam(Cplex.Param.MIP.Pool.Replace, 2); // 解替换策略
                                                                      //设置为0时，新找到的解将替换解池中的解，只保留最优解。
                                                                      //设置为1时，新找到的解将替换解池中的解，但保留一些次优解。
                                                                      //设置为2时，新找到的解将与解池中的解进行比较，如果新解更好，则替换解池中的解。

            ssmaster.model.SetParam(Cplex.Param.MIP.Pool.Intensity, 0); // 解池强度 
                                                                        //较大的值表示解池中保留的解较多，可以存储更多的解。
                                                                        //较小的值表示解池中保留的解较少，仅保留最好的几个解。
                                                                        // 设置变量选择策略为full-strong branching            

            ssmaster.model.Use(new SCNRLazyConsCallback (ssmaster, data, y_sol));

            ssmaster.model.Solve();
            Cplex.Status modelsts = ssmaster.model.GetStatus();
            scenarioCost = ssmaster.model.GetObjValue();
            double bestScenarioCost = ssmaster.model.GetBestObjValue();

            Console.WriteLine($"Optimal objective cost: {scenarioCost}");
            Cplex.Status status = ssmaster.model.GetStatus();
            int[] u_sol = new int[data.sectionSize];
            for (int a = 0; a < data.sectionSize; a++)
            {
                double val = ssmaster.model.GetValue(ssmaster.u_a[a]);
                if (val > 0.5)
                {
                    u_sol[a] = 1;
                }
            }
            //in case GetObjValue()!= GetBestObjValue()
            if (bestScenarioCost > scenarioCost)
            {
                scenarioCost = ssmaster.bestObjVal;
                xI_sols_pool.Add(u_sol);
                u_sol = ssmaster.bestFeasibleSolution;
            }

            int numSols = ssmaster.model.GetSolnPoolNsolns();

            if (numSols >= 1)
            {
                List<int[]> candidatesolset = new List<int[]>();
                List<double> costset = new List<double>();

                for (int n = 0; n < numSols; n++)
                {
                    bool opttest = true;
                    int[] feasiblesolution = new int[data.sectionSize];
                    for (int a = 0; a < data.sectionSize; a++)
                    {
                        if (ssmaster.model.GetValue(ssmaster.u_a[a], n) > 0.5)
                        {
                            feasiblesolution[a] = 1;
                        }
                        if (feasiblesolution[a] != u_sol[a])
                            opttest = false;
                    }
                    double objval = ssmaster.model.GetValue(ssmaster.eta, n);

                    if (!opttest)
                    {
                        costset.Add(objval);
                        candidatesolset.Add(feasiblesolution);
                    }
                }

                int counter = 0;
                while (true)
                {
                    if (counter >= data.maxMultiCuts) break;

                    if (costset.Count < 1) break;


                    int ind = costset.IndexOf(costset.Max());
                    xI_sols_pool.Add(candidatesolset[ind]);

                    candidatesolset.RemoveAt(ind); costset.RemoveAt(ind);
                    counter++;
                }
            }

            return ssmaster;
            #endregion
        }
    }
    class CandCG
    {
        Data data = new Data();
        Solution solution = new Solution();
        CCGSub subproblem;
        CCGMaster masterproblem = new CCGMaster();

        int iter_l = 0;

        public void ParallelCCG()
        {
            #region solution process
            //generate master problem and subproblem
            List<Path> feedbacktomaster = new List<Path>();
            double eta_value = 0, mastercost = 0, UB = float.MaxValue, LB = float.MinValue;
            bool isparrallel = true;
            int cutting_flag = 0;// add the cutting plane
            int[] y_solution = new int[data.DCSize];
            
            masterproblem.InitializeCCGMaster();
            masterproblem.model.SetParam(Cplex.Param.Threads, 1);

            int number_of_opens = 0;
            for (int j = 0; j < data.DCSize; j++)
            {
                y_solution[j] = 1;
                number_of_opens++;
            }

            bool feasible_flag = false;
            double worst_scenario_cost = float.MinValue;
            int[] worst_u_a = new int[data.sectionSize];
            int[] section_to_link_state = new int[data.linkSize];

            List<int[]> scenarios = new List<int[]>();
            List<Cplex> subproblemlist = new List<Cplex>();

            int number_of_scenarios = 0;
            // Find the worst scenario in the context of the current solution, i.e., two cases
            //case 1: the number of open DCents less than or equal to the predefined threshold
            worst_scenario_cost = float.MinValue;
            feasible_flag = true;//indicate whether the subproblem is feasible or not
            subproblemlist = new List<Cplex>();

            List<int> open_DCent_No_list = new List<int>();
            for (int j = 0; j < data.sectionSize; j++)
            {
                open_DCent_No_list.Add(j);
            }

            // All scenarios
            List<List<int>> all_scenarios = new List<List<int>>();

            //no decents is destroyed                
            //enumarate all scenarios
            for (int i = 1; i <= data.max_dstroyed_DCs; i++)
            {
                List<List<int>> temp_scenarios = Combination.Combine(open_DCent_No_list, i);

                for (int l = 0; l < temp_scenarios.Count; l++)
                {
                    List<int> tem_destroy_solution = temp_scenarios[l].ToList();


                    all_scenarios.Add(tem_destroy_solution.ToList());
                }
            }

            number_of_scenarios = all_scenarios.Count;

            for (int l = 0; l < all_scenarios.Count; l++)
            {
                List<int> temp_solution = all_scenarios[l].ToList();

                int[] Dcent_state = new int[data.sectionSize];

                for (int j = 0; j < data.sectionSize; j++)
                {
                    Dcent_state[j] = 0;//initialization;                        
                }

                for (int j = 0; j < temp_solution.Count; j++)
                {
                    int dcentno = temp_solution[j];
                    Dcent_state[dcentno] = 1;//diruption                       
                }
                for (int s = 0; s < data.sectionSize; s++)
                {
                    //transfer section state to link state
                    if (Dcent_state[s] == 1)
                    {
                        Section section = data.sectionlist[s];
                        for (int a = 0; a < section.CorrespondLinks.Count; a++)
                        {
                            section_to_link_state[section.CorrespondLinks[a].LinkSeqNo] = 1;
                        }
                    }
                }
                subproblem.GenCCGSubproblem(y_solution, section_to_link_state);
                //subproblem.model.ExportModel("CCGsubproblem.lp");
                //collect subproblem                  
                subproblemlist.Add(subproblem.model);
                //collect scenarios
                scenarios.Add(Dcent_state);
                if (!isparrallel)
                {
                    if (subproblem.model.Solve())
                    {
                        if (subproblem.model.GetObjValue() > worst_scenario_cost)
                        {
                            worst_scenario_cost = subproblem.model.GetObjValue();
                            worst_u_a = Dcent_state.ToArray();
                        }
                    }
                    else
                    {
                        feasible_flag = false;
                        worst_u_a = Dcent_state.ToArray();
                        break;
                    }
                }
            }

            if (isparrallel)
            {
                bool infeasibleFound = false;  // 共享标志
                object lockObj = new object();  // 锁对象，用于线程同步

                var options = new ParallelOptions() { MaxDegreeOfParallelism = Environment.ProcessorCount };
                // 使用Parallel.ForEach来并行处理
                Parallel.For(0, subproblemlist.Count, options, (i, state) =>
                {
                    if (infeasibleFound)
                    {
                        state.Stop();
                        return;
                    }
                    subproblemlist[i].Solve();
                    Cplex.Status isfeasible = subproblemlist[i].GetStatus();
                    if (isfeasible == Cplex.Status.Infeasible)
                    {
                        lock (lockObj)
                        {
                            if (!infeasibleFound)  // 再次检查以避免竞态条件
                            {
                                infeasibleFound = true;
                                feasible_flag = false;
                                worst_u_a = scenarios[i].ToArray();
                            }
                        }
                        state.Stop();
                    }
                    else
                    {
                        double objValue = subproblemlist[i].GetObjValue();
                        lock (lockObj)
                        {
                            if (objValue > worst_scenario_cost)
                            {
                                worst_scenario_cost = objValue;
                                worst_u_a = scenarios[i].ToArray();
                            }
                        }
                    }
                });
            }

            if (feasible_flag)
            {
                //find the worst scenario                    
                UB = Math.Min(UB, worst_scenario_cost + mastercost - eta_value);
            }
            if (Program.CCG_debug == 1)
            {
                Program.CCG_debugfile.Write("-----------------------------------------------------------------------------------------------------------" +
                    "--------------------------------------------------------------------------------------------------------------\r\n");
                Program.CCG_debugfile.WriteLine("Iter:{0}, upper bound: {1}, lower bound: {2}, gap: {3} value of eta: {4}, master value: {5}, subproblem value: {6}, " +
                    "number of scenarios: {7}, feasible state: {8}",
                    iter_l, UB, LB, (UB - LB) / UB * 100, eta_value, mastercost, worst_scenario_cost, number_of_scenarios, feasible_flag);
                Program.CCG_debugfile.Write("-----------------------------------------------------------------------------------------------------------" +
                    "--------------------------------------------------------------------------------------------------------------\r\n");
            }
            double eps = 0.001;
            iter_l = -1;
            while (UB - LB > eps)
            {
                iter_l++;
                if (feasible_flag)// existing a worst scenario
                {
                    section_to_link_state = new int[data.linkSize];
                    cutting_flag = 1;//add the cutting plane
                    for (int s = 0; s < data.sectionSize; s++)
                    {
                        //transfer section state to link state
                        if (worst_u_a[s] == 1)
                        {
                            Section section = data.sectionlist[s];
                            for (int a = 0; a < section.CorrespondLinks.Count; a++)
                            {
                                section_to_link_state[section.CorrespondLinks[a].LinkSeqNo] = 1;
                            }
                        }
                    }
                    
                    masterproblem.u_l_a_to_link_state.Add(section_to_link_state);

                    masterproblem.GenCCGMasterproblem(cutting_flag);// master(cutting_flag, paths_at_iterl);
                    masterproblem.model.ExportModel("CCGMaster.lp");
                    masterproblem.model.Solve();
                    LB = masterproblem.model.ObjValue;
                    mastercost = masterproblem.model.ObjValue;
                    eta_value = masterproblem.model.GetValue(masterproblem.omega);

                    //update the facility location decision
                    y_solution = new int[data.DCSize];

                    number_of_opens = 0;
                    for (int j = 0; j < data.DCSize; j++)
                    {
                        y_solution[j] = 0;
                        double y_j_value = masterproblem.model.GetValue(masterproblem.y_j[j]);

                        if (y_j_value > 0.5)
                        {
                            y_solution[j] = 1;
                            number_of_opens++;

                        }
                    }
                    
                    feasible_flag = false;

                    worst_u_a = new int[data.sectionSize];
                    section_to_link_state = new int[data.linkSize];

                    scenarios = new List<int[]>();

                    // Find the worst scenario in the context of the current solution, i.e., two cases
                    //case 1: the number of open DCents less than or equal to the predefined threshold
                    worst_scenario_cost = float.MinValue;
                    feasible_flag = true;//indicate whether the subproblem is feasible or not
                    subproblemlist = new List<Cplex>();

                    open_DCent_No_list = new List<int>();
                    for (int j = 0; j < data.sectionSize; j++)
                    {
                        open_DCent_No_list.Add(j);
                    }

                    // All scenarios
                    all_scenarios = new List<List<int>>();

                    //no decents is destroyed                
                    //enumarate all scenarios
                    for (int i = 1; i <= data.max_dstroyed_DCs; i++)
                    {
                        List<List<int>> temp_scenarios = Combination.Combine(open_DCent_No_list, i);

                        for (int l = 0; l < temp_scenarios.Count; l++)
                        {
                            List<int> tem_destroy_solution = temp_scenarios[l].ToList();


                            all_scenarios.Add(tem_destroy_solution.ToList());
                        }
                    }

                    number_of_scenarios = all_scenarios.Count;

                    for (int l = 0; l < all_scenarios.Count; l++)
                    {
                        List<int> temp_solution = all_scenarios[l].ToList();

                        int[] Dcent_state = new int[data.sectionSize];

                        for (int j = 0; j < data.sectionSize; j++)
                        {
                            Dcent_state[j] = 0;//initialization;                        
                        }

                        for (int j = 0; j < temp_solution.Count; j++)
                        {
                            int dcentno = temp_solution[j];
                            Dcent_state[dcentno] = 1;//diruption                       
                        }
                        for (int s = 0; s < data.sectionSize; s++)
                        {
                            //transfer section state to link state
                            if (Dcent_state[s] == 1)
                            {
                                Section section = data.sectionlist[s];
                                for (int a = 0; a < section.CorrespondLinks.Count; a++)
                                {
                                    section_to_link_state[section.CorrespondLinks[a].LinkSeqNo] = 1;
                                }
                            }
                        }
                        subproblem.GenCCGSubproblem(y_solution, section_to_link_state);
                        //subproblem.model.ExportModel("CCGsubproblem.lp");
                        //collect subproblem                  
                        subproblemlist.Add(subproblem.model);
                        //collect scenarios
                        scenarios.Add(Dcent_state);
                        if (!isparrallel)
                        {
                            if (subproblem.model.Solve())
                            {
                                if (subproblem.model.GetObjValue() > worst_scenario_cost)
                                {
                                    worst_scenario_cost = subproblem.model.GetObjValue();
                                    worst_u_a = Dcent_state.ToArray();
                                }
                            }
                            else
                            {
                                feasible_flag = false;
                                worst_u_a = Dcent_state.ToArray();
                                break;
                            }
                        }
                    }

                    if (isparrallel)
                    {
                        bool infeasibleFound = false;  // 共享标志
                        object lockObj = new object();  // 锁对象，用于线程同步

                        var options = new ParallelOptions() { MaxDegreeOfParallelism = Environment.ProcessorCount };
                        // 使用Parallel.ForEach来并行处理
                        Parallel.For(0, subproblemlist.Count, options, (i, state) =>
                        {
                            if (infeasibleFound)
                            {
                                state.Stop();
                                return;
                            }
                            subproblemlist[i].Solve();
                            Cplex.Status isfeasible = subproblemlist[i].GetStatus();
                            if (isfeasible == Cplex.Status.Infeasible)
                            {
                                lock (lockObj)
                                {
                                    if (!infeasibleFound)  // 再次检查以避免竞态条件
                                    {
                                        infeasibleFound = true;
                                        feasible_flag = false;
                                        worst_u_a = scenarios[i].ToArray();
                                    }
                                }
                                state.Stop();
                            }
                            else
                            {
                                double objValue = subproblemlist[i].GetObjValue();
                                lock (lockObj)
                                {
                                    if (objValue > worst_scenario_cost)
                                    {
                                        worst_scenario_cost = objValue;
                                        worst_u_a = scenarios[i].ToArray();
                                    }
                                }
                            }
                        });
                    }
                    //find the worst scenario
                    if (feasible_flag == true)
                    {
                        //find the worst scenario
                        UB = Math.Min(UB, worst_scenario_cost + mastercost - eta_value);
                    }
                    if (Program.CCG_debug == 1)
                    {
                        Program.CCG_debugfile.Write("-----------------------------------------------------------------------------------------------------------" +
                            "--------------------------------------------------------------------------------------------------------------\r\n");
                        Program.CCG_debugfile.WriteLine("Iter:{0}, upper bound: {1}, lower bound: {2}, gap: {3} value of eta: {4}, master value: {5}, subproblem value: {6}, " +
                            "number of scenarios: {7}, feasible state: {8}",
                            iter_l, UB, LB, (UB - LB) / UB * 100, eta_value, mastercost, worst_scenario_cost, number_of_scenarios, feasible_flag);
                        Program.CCG_debugfile.Write("-----------------------------------------------------------------------------------------------------------" +
                            "--------------------------------------------------------------------------------------------------------------\r\n");
                    }
                    
                }
                else //unbounded
                {
                    cutting_flag = 0;//do not add the cutting plane
                    for (int s = 0; s < data.sectionSize; s++)
                    {
                        //transfer section state to link state
                        if (worst_u_a[s] == 1)
                        {
                            Section section = data.sectionlist[s];
                            for (int a = 0; a < section.CorrespondLinks.Count; a++)
                            {
                                section_to_link_state[section.CorrespondLinks[a].LinkSeqNo] = 1;
                            }
                        }
                    }
                    masterproblem.u_l_a_to_link_state.Add(section_to_link_state);

                    masterproblem.GenCCGMasterproblem(cutting_flag);
                    masterproblem.model.ExportModel("CCGMaster.lp");
                    masterproblem.model.Solve();

                    LB = masterproblem.model.ObjValue;
                    mastercost = masterproblem.model.ObjValue;
                    eta_value = masterproblem.model.GetValue(masterproblem.omega);

                    //update the facility location decision
                    y_solution = new int[data.DCSize];

                    number_of_opens = 0;
                    for (int j = 0; j < data.DCSize; j++)
                    {
                        y_solution[j] = 0;
                        double y_j_value = masterproblem.model.GetValue(masterproblem.y_j[j]);

                        if (y_j_value > 0.5)
                        {
                            y_solution[j] = 1;
                            number_of_opens++;

                        }

                    }
                    
                    feasible_flag = false;

                    worst_u_a = new int[data.sectionSize];
                    section_to_link_state = new int[data.linkSize];
                    scenarios = new List<int[]>();

                    // Find the worst scenario in the context of the current solution, i.e., two cases
                    //case 1: the number of open DCents less than or equal to the predefined threshold
                    worst_scenario_cost = float.MinValue;
                    feasible_flag = true;//indicate whether the subproblem is feasible or not
                    subproblemlist = new List<Cplex>();

                    open_DCent_No_list = new List<int>();
                    for (int j = 0; j < data.sectionSize; j++)
                    {
                        open_DCent_No_list.Add(j);
                    }

                    // All scenarios
                    all_scenarios = new List<List<int>>();

                    //no decents is destroyed                
                    //enumarate all scenarios
                    for (int i = 1; i <= data.max_dstroyed_DCs; i++)
                    {
                        List<List<int>> temp_scenarios = Combination.Combine(open_DCent_No_list, i);

                        for (int l = 0; l < temp_scenarios.Count; l++)
                        {
                            List<int> tem_destroy_solution = temp_scenarios[l].ToList();

                            all_scenarios.Add(tem_destroy_solution.ToList());
                        }
                    }

                    number_of_scenarios = all_scenarios.Count;

                    for (int l = 0; l < all_scenarios.Count; l++)
                    {
                        List<int> temp_solution = all_scenarios[l].ToList();

                        int[] Dcent_state = new int[data.sectionSize];

                        for (int j = 0; j < data.sectionSize; j++)
                        {
                            Dcent_state[j] = 0;//initialization;                        
                        }

                        for (int j = 0; j < temp_solution.Count; j++)
                        {
                            int dcentno = temp_solution[j];
                            Dcent_state[dcentno] = 1;//diruption                       
                        }
                        for (int s = 0; s < data.sectionSize; s++)
                        {
                            //transfer section state to link state
                            if (Dcent_state[s] == 1)
                            {
                                Section section = data.sectionlist[s];
                                for (int a = 0; a < section.CorrespondLinks.Count; a++)
                                {
                                    section_to_link_state[section.CorrespondLinks[a].LinkSeqNo] = 1;
                                }
                            }
                        }
                        subproblem.GenCCGSubproblem(y_solution, section_to_link_state);
                        //subproblem.model.ExportModel("CCGsubproblem.lp");
                        //collect subproblem                  
                        subproblemlist.Add(subproblem.model);
                        //collect scenarios
                        scenarios.Add(Dcent_state);
                        if (!isparrallel)
                        {
                            if (subproblem.model.Solve())
                            {
                                if (subproblem.model.GetObjValue() > worst_scenario_cost)
                                {
                                    worst_scenario_cost = subproblem.model.GetObjValue();
                                    worst_u_a = Dcent_state.ToArray();
                                }
                            }
                            else
                            {
                                feasible_flag = false;
                                worst_u_a = Dcent_state.ToArray();
                                break;
                            }
                        }
                    }

                    if (isparrallel)
                    {
                        bool infeasibleFound = false;  // 共享标志
                        object lockObj = new object();  // 锁对象，用于线程同步

                        var options = new ParallelOptions() { MaxDegreeOfParallelism = Environment.ProcessorCount };
                        // 使用Parallel.ForEach来并行处理
                        Parallel.For(0, subproblemlist.Count, options, (i, state) =>
                        {
                            if (infeasibleFound)
                            {
                                state.Stop();
                                return;
                            }
                            subproblemlist[i].Solve();
                            Cplex.Status isfeasible = subproblemlist[i].GetStatus();
                            if (isfeasible == Cplex.Status.Infeasible)
                            {
                                lock (lockObj)
                                {
                                    if (!infeasibleFound)  // 再次检查以避免竞态条件
                                    {
                                        infeasibleFound = true;
                                        feasible_flag = false;
                                        worst_u_a = scenarios[i].ToArray();
                                    }
                                }
                                state.Stop();
                            }
                            else
                            {
                                double objValue = subproblemlist[i].GetObjValue();
                                lock (lockObj)
                                {
                                    if (objValue > worst_scenario_cost)
                                    {
                                        worst_scenario_cost = objValue;
                                        worst_u_a = scenarios[i].ToArray();
                                    }
                                }
                            }
                        });
                    }
                    //find the worst scenario
                    if (feasible_flag == true)
                    {
                        //find the worst scenario
                        UB = Math.Min(UB, worst_scenario_cost + mastercost - eta_value);
                    }
                    if (Program.CCG_debug == 1)
                    {
                        Program.CCG_debugfile.Write("-----------------------------------------------------------------------------------------------------------" +
                            "--------------------------------------------------------------------------------------------------------------\r\n");
                        Program.CCG_debugfile.WriteLine("Iter:{0}, upper bound: {1}, lower bound: {2}, gap: {3} value of eta: {4}, master value: {5}, subproblem value: {6}, " +
                            "number of scenarios: {7}, feasible state: {8}",
                            iter_l, UB, LB, (UB - LB) / UB * 100, eta_value, mastercost, worst_scenario_cost, number_of_scenarios, feasible_flag);
                        Program.CCG_debugfile.Write("-----------------------------------------------------------------------------------------------------------" +
                            "--------------------------------------------------------------------------------------------------------------\r\n");
                    }
                    continue;
                }
            }

            //solution.write_solution(masterproblem);
            Console.WriteLine("UB:{0}, LB:{1}", UB, LB);
            #endregion
        }
        public void CCGwBAC()
        {
            #region solution process
            //generate master problem and subproblem
            List<Path> feedbacktomaster = new List<Path>();
            double eta_value = 0, mastercost = 0, UB = float.MaxValue, LB = float.MinValue;
            
            int cutting_flag = 0;// add the cutting plane
            int[] y_solution = new int[data.DCSize];

            masterproblem.InitializeCCGMaster();
            masterproblem.model.SetParam(Cplex.Param.Threads, 1);

            int number_of_opens = 0;
            for (int j = 0; j < data.DCSize; j++)
            {
                y_solution[j] = 1;
                if (number_of_opens == 2)
                    break;
                number_of_opens++;
            }

            bool feasible_flag = false;
            double worst_scenario_cost = float.MinValue;
            
            int[] worst_u_a = new int[data.sectionSize];
            int[] section_to_link_state = new int[data.linkSize];

            CACGBD colBDGen = new CACGBD();

            SCNRMaster submaster = colBDGen.BranchandCutForWorstCsenario(y_solution);

            Cplex.Status cur_status = submaster.model.GetStatus();
            Console.WriteLine(submaster.model.GetObjValue());
            if (cur_status == Cplex.Status.Optimal || cur_status == Cplex.Status.Feasible)
            {
                feasible_flag = true;
                worst_scenario_cost = submaster.model.GetObjValue();
                for (int l = 0; l < data.sectionSize; l++)
                {
                    var u_val = submaster.model.GetValue(submaster.u_a[l]);

                    if(u_val > 0.5)
                    {
                        worst_u_a[l] = 1;
                    }
                    
                }
            }
            else if (cur_status == Cplex.Status.Unbounded)
            {
                feasible_flag = false;
                
                for (int l = 0; l < data.sectionSize; l++)
                {
                    var u_val = submaster.model.GetValue(submaster.u_a[l]);

                    if (u_val > 0.5)
                    {
                        worst_u_a[l] = 1;
                    }
                }
            }

            if (Program.CCG_debug == 1)
            {
                Program.CCG_debugfile.Write("-----------------------------------------------------------------------------------------------------------" +
                    "--------------------------------------------------------------------------------------------------------------\r\n");
                Program.CCG_debugfile.WriteLine("Iter:{0}, upper bound: {1}, lower bound: {2}, gap: {3} value of eta: {4}, master value: {5}, subproblem value: {6}, " +
                    "feasible state: {7}",
                    iter_l, UB, LB, (UB - LB) / UB * 100, eta_value, mastercost, worst_scenario_cost, feasible_flag);
                Program.CCG_debugfile.Write("-----------------------------------------------------------------------------------------------------------" +
                    "--------------------------------------------------------------------------------------------------------------\r\n");
            }

            double eps = 0.001;
            iter_l = -1;
            while (UB - LB > eps)
            {
                iter_l++;
                if (feasible_flag)// existing a worst scenario
                {
                    cutting_flag = 1;//add the cutting plane
                    section_to_link_state = new int[data.linkSize];

                    for (int s = 0; s < data.sectionSize; s++)
                    {
                        //transfer section state to link state
                        if(worst_u_a[s] == 1)
                        {
                            Section section = data.sectionlist[s];
                            for (int l = 0; l < section.CorrespondLinks.Count; l++)
                            {
                                section_to_link_state[section.CorrespondLinks[l].LinkSeqNo] = 1;
                            }
                        }                        
                    }
                    masterproblem.u_l_a_to_link_state.Add(section_to_link_state);

                    masterproblem.GenCCGMasterproblem(cutting_flag);// master(cutting_flag, paths_at_iterl);
                    masterproblem.model.ExportModel("CCGMaster.lp");
                    masterproblem.model.Solve();
                    LB = masterproblem.model.ObjValue;
                    mastercost = masterproblem.model.ObjValue;
                    eta_value = masterproblem.model.GetValue(masterproblem.omega);

                    //update the facility location decision
                    y_solution = new int[data.DCSize];

                    number_of_opens = 0;
                    for (int j = 0; j < data.DCSize; j++)
                    {
                        y_solution[j] = 0;
                        double y_j_value = masterproblem.model.GetValue(masterproblem.y_j[j]);

                        if (y_j_value > 0.5)
                        {
                            y_solution[j] = 1;
                            number_of_opens++;

                        }
                    }
                    for (int l = 0; l < masterproblem.x_l_r.Count; l++)
                    {
                        for (int r = 0; r < data.pathSize; r++)
                        {
                            double x_r_value = masterproblem.model.GetValue(masterproblem.x_l_r[l][r]);
                            if (x_r_value > 0)
                            {
                                if (data.pathlist[r].Contains(44))
                                    Console.WriteLine();
                                Console.WriteLine($"iteration:{l},xrvalue:{x_r_value},pathnumber:{r},startnode:{Program.g_DCent_list[data.s_r[r]].DCenterID},endnode:{Program.g_node_list[data.e_r[r]].NodeID}");
                            }
                           
                        }
                        Console.WriteLine();
                    }
                    for (int l = 0; l < masterproblem.z_l_j.Count; l++)
                    {
                        for (int j = 0; j < data.DCSize; j++)
                        {
                            double z_j_value = masterproblem.model.GetValue(masterproblem.z_l_j[l][j]);
                            if (z_j_value > 0)
                            {
                                Console.WriteLine($"iteration:{l},jvalue:{j},zjvalue:{z_j_value}");
                            }
                        }
                        Console.WriteLine();
                    }

                    feasible_flag = false;
                    
                    worst_u_a = new int[data.sectionSize];
                    
                    // Find the worst scenario in the context of the current solution

                    submaster = colBDGen.BranchandCutForWorstCsenario(y_solution);
                    cur_status = submaster.model.GetStatus();
                    Console.WriteLine(submaster.model.GetObjValue());

                    if (cur_status == Cplex.Status.Optimal || cur_status == Cplex.Status.Feasible)
                    {
                        feasible_flag = true;
                        worst_scenario_cost = submaster.model.GetObjValue();
                        for (int l = 0; l < data.sectionSize; l++)
                        {
                            var u_val = submaster.model.GetValue(submaster.u_a[l]);

                            if (u_val > 0.5)
                            {
                                worst_u_a[l] = 1;
                            }
                        }
                    }
                    else if (cur_status == Cplex.Status.Unbounded)
                    {
                        feasible_flag = false;

                        for (int l = 0; l < data.sectionSize; l++)
                        {
                            var u_val = submaster.model.GetValue(submaster.u_a[l]);

                            if (u_val > 0.5)
                            {
                                worst_u_a[l] = 1;
                            }
                        }
                    }
                    if (feasible_flag)
                    {
                        //find the worst scenario                    
                        UB = Math.Min(UB, worst_scenario_cost + mastercost - eta_value);
                    }
                    for (int r = 0; r < data.pathSize; r++)
                    {
                        if (data.pathlist[r].Contains(111))
                            Console.WriteLine($"{r},{data.s_r[r]},{data.e_r[r]}");
                    }

                    if (Program.CCG_debug == 1)
                    {
                        Program.CCG_debugfile.Write("-----------------------------------------------------------------------------------------------------------" +
                            "--------------------------------------------------------------------------------------------------------------\r\n");
                        Program.CCG_debugfile.WriteLine("Iter:{0}, upper bound: {1}, lower bound: {2}, gap: {3} value of eta: {4}, master value: {5}, subproblem value: {6}, " +
                            "feasible state: {7}",
                            iter_l, UB, LB, (UB - LB) / UB * 100, eta_value, mastercost, worst_scenario_cost, feasible_flag);
                        Program.CCG_debugfile.Write("-----------------------------------------------------------------------------------------------------------" +
                            "--------------------------------------------------------------------------------------------------------------\r\n");
                    }
                }
                else //unbounded
                {
                    cutting_flag = 0;//do not add the cutting plane

                    for (int s = 0; s < data.sectionSize; s++)
                    {
                        //transfer section state to link state
                        if (worst_u_a[s] == 1)
                        {
                            Section section = data.sectionlist[s];
                            for (int l = 0; l < section.CorrespondLinks.Count; l++)
                            {
                                section_to_link_state[section.CorrespondLinks[l].LinkSeqNo] = 1;
                            }
                        }
                    }
                    masterproblem.u_l_a_to_link_state.Add(section_to_link_state);

                    masterproblem.GenCCGMasterproblem(cutting_flag);
                    masterproblem.model.ExportModel("CCGMaster.lp");
                    masterproblem.model.Solve();

                    LB = masterproblem.model.ObjValue;
                    mastercost = masterproblem.model.ObjValue;
                    eta_value = masterproblem.model.GetValue(masterproblem.omega);

                    //update the facility location decision
                    y_solution = new int[data.DCSize];

                    number_of_opens = 0;
                    for (int j = 0; j < data.DCSize; j++)
                    {
                        y_solution[j] = 0;
                        double y_j_value = masterproblem.model.GetValue(masterproblem.y_j[j]);

                        if (y_j_value > 0.5)
                        {
                            y_solution[j] = 1;
                            number_of_opens++;

                        }

                    }
                    feasible_flag = false;

                    worst_u_a = new int[data.sectionSize];
                    
                    // Find the worst scenario in the context of the current solution

                    submaster = colBDGen.BranchandCutForWorstCsenario(y_solution);
                    cur_status = submaster.model.GetStatus();

                    if (cur_status == Cplex.Status.Optimal || cur_status == Cplex.Status.Feasible)
                    {
                        feasible_flag = true;
                        worst_scenario_cost = submaster.model.GetObjValue();
                        for (int l = 0; l < data.sectionSize; l++)
                        {
                            var u_val = submaster.model.GetValue(submaster.u_a[l]);

                            if (u_val > 0.5)
                            {
                                worst_u_a[l] = 1;
                            }
                        }
                    }
                    else if (cur_status == Cplex.Status.Unbounded)
                    {
                        feasible_flag = false;

                        for (int l = 0; l < data.sectionSize; l++)
                        {
                            var u_val = submaster.model.GetValue(submaster.u_a[l]);

                            if (u_val > 0.5)
                            {
                                worst_u_a[l] = 1;
                            }
                        }
                    }
                    if (feasible_flag)
                    {
                        //find the worst scenario                    
                        UB = Math.Min(UB, worst_scenario_cost + mastercost - eta_value);
                    }
                    if (Program.CCG_debug == 1)
                    {
                        Program.CCG_debugfile.Write("-----------------------------------------------------------------------------------------------------------" +
                            "--------------------------------------------------------------------------------------------------------------\r\n");
                        Program.CCG_debugfile.WriteLine("Iter:{0}, upper bound: {1}, lower bound: {2}, gap: {3} value of eta: {4}, master value: {5}, subproblem value: {6}, " +
                            "feasible state: {7}",
                            iter_l, UB, LB, (UB - LB) / UB * 100, eta_value, mastercost, worst_scenario_cost, feasible_flag);
                        Program.CCG_debugfile.Write("-----------------------------------------------------------------------------------------------------------" +
                            "--------------------------------------------------------------------------------------------------------------\r\n");
                    }
                }
            }

            //solution.write_solution(masterproblem);
            Console.WriteLine("UB:{0}, LB:{1}", UB, LB);
            #endregion
        }        
    }
}
