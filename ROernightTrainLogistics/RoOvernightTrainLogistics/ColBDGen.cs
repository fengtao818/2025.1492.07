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
using System.Numerics;
using System.Reflection;

namespace RoOvernightTrainLogistics
{
    class DataStructure
    {
        #region
        public double righthandsideU;

        public double[] x_r_val;        
        public double[] v_i_val;
        public double[] delta_r_val;
        public double[] w_j_val;
        public double[] pi_a_val;

        public double[,] MXi_val;
        public DataStructure(Data data)
        {
            x_r_val = new double[data.pathSize];
            v_i_val = new double[data.nodeSize];
            MXi_val = new double[data.nodeSize, data.DCSize];
            delta_r_val = new double[data.pathSize];
            w_j_val = new double[data.DCSize];
            pi_a_val = new double[data.linkSize];
        }
        #endregion
    }
    internal class BendersLazyConsCallback : Cplex.LazyConstraintCallback
    {
        internal BendersMaster Bendmaster;
        internal Data data;
        
        internal int count;
        
        internal CCGSub ccgsub;
        internal DualCCGSub dccgsub;
        internal BendersDualSub bddsub;
        Random rand;
        bool updatePareto;
        List<int[]> scenariolist;
        
        public BendersLazyConsCallback(BendersMaster Bendmaster, Data dt)
        {
            this.Bendmaster = Bendmaster;
            
            data = dt;
            count = 0;
            
            ccgsub = new CCGSub();
            ccgsub.GenCCGSubproblem(new int[data.DCSize], new int[data.DCSize]);
            ccgsub.model.SetParam(Cplex.Param.Preprocessing.Presolve, false);            
            

            bddsub = new BendersDualSub();
            bddsub.GendeltaUBsubDual(new int[data.DCSize], new int[data.DCSize]);
            bddsub.model.SetParam(Cplex.Param.Preprocessing.Presolve, false);

            dccgsub = new DualCCGSub();
            dccgsub.GenDualofCCGsub(new int[data.DCSize], new double[data.pathSize]);
            dccgsub.model.SetParam(Cplex.Param.MIP.Pool.Capacity, data.maxMultiCuts); 
            dccgsub.model.SetParam(Cplex.Param.MIP.Pool.Replace, 2);
            dccgsub.model.SetParam(Cplex.Param.MIP.Pool.Intensity, 0);

            updatePareto = true;

            scenariolist = new List<int[]>();

            rand = new Random();

        }
        public override void Main()
        {
            #region BranchandCheck
            
            int[] y_solution = new int[data.DCSize];
            int[] worst_u_j = new int[data.DCSize];           
            List<int[]> multicutslist = new List<int[]>();

            Stopwatch subproblemtime = new Stopwatch();
            Stopwatch multicutstime = new Stopwatch();
            Stopwatch paretotime = new Stopwatch();
            Stopwatch findw0rstscenario = new Stopwatch();

            double objectivecost = GetObjValue();
            double omega = GetValue(Bendmaster.omega);
            double subproblemtimerecord = 0;
            double multicutstimerecord = 0;
            double paretotimerecord = 0;
            double findw0rstscenariorecord = 0;
            
            int numberofopens = 0;
            for (int j = 0; j < data.DCSize; j++)
            {
                
                var yval = GetValue(Bendmaster.y_j[j]);
                
                if (yval > 0.5)
                {
                    y_solution[j] = 1;
                    numberofopens++;
                }                
            }
            
            count++;
            
            CCGSub ccgsubproblem = new CCGSub();
            List<Cplex> subproblemlist = new List<Cplex>();
            List<int[]> scenarios = new List<int[]>();
            List<double> costlist = new List<double>();

            object lockObject = new object();
            double worst_scenario_cost = float.MinValue;
            
            
            bool issubprofeasi = true;
            if (data.BACforworstscenario)
            {
                findw0rstscenario.Start();
                Dictionary<int, double> gvval = new Dictionary<int, double>();
                CACGBD colBDGen = new CACGBD();
                SCNRMaster submaster = colBDGen.BranchandCutForWorstCsenario(y_solution);

                
                issubprofeasi = submaster.feasiblestatus;
                
                if (!issubprofeasi)
                {
                    worst_u_j = submaster.bestFeasibleSolution.ToArray();

                }
                else
                {
                    worst_u_j = new int[data.DCSize];
                    for (int l = 0; l < data.DCSize; l++)
                    {
                        var u_val = submaster.model.GetValue(submaster.u_j[l]);
                        if (u_val > 0.5)
                        {
                            worst_u_j[l] = 1;
                        }
                    }
                    multicutslist = colBDGen.xI_sols_pool.ToList();
                    
                    ccgsub.ResetOBJaRHS_XI(y_solution, worst_u_j);
                    ccgsub.model.Solve();
                    worst_scenario_cost = ccgsub.model.GetObjValue();
                }
                submaster.model.End();
                findw0rstscenario.Stop(); findw0rstscenariorecord = findw0rstscenario.ElapsedMilliseconds;               
            }
            else
            {
                findw0rstscenario.Start();

                
                
                
                bddsub.ResetBDDSubObjforDelta(y_solution, y_solution);
                bddsub.model.Solve();

                double[] deltaval = new double[data.pathSize];
                deltaval = bddsub.model.GetValues(bddsub.delta_r);

                
                dccgsub.ResetdualccgSubObjcons(y_solution, deltaval);
                dccgsub.model.Solve();

                worst_u_j = new int[data.DCSize];
                for (int l = 0; l < data.DCSize; l++)
                {
                    var u_val = dccgsub.model.GetValue(dccgsub.u_j[l]);

                    if (u_val > 0.5)
                    {
                        worst_u_j[l] = 1;
                    }                    
                }
                issubprofeasi = true;
                
                ccgsub.ResetOBJaRHS_XI(y_solution, worst_u_j);
                issubprofeasi = ccgsub.model.Solve();
                worst_scenario_cost = ccgsub.model.GetObjValue();

                int numSols = dccgsub.model.GetSolnPoolNsolns();

                if (numSols >= 1)
                {
                    List<int[]> candidatesolset = new List<int[]>();
                    List<double> costset = new List<double>();

                    for (int n = 0; n < numSols; n++)
                    {
                        bool opttest = true;
                        int[] feasiblesolution = new int[data.DCSize];
                        for (int j = 0; j < data.DCSize; j++)
                        {
                            if (dccgsub.model.GetValue(dccgsub.u_j[j], n) > 0.5)
                            {
                                feasiblesolution[j] = 1;
                            }
                            if (feasiblesolution[j] != worst_u_j[j])
                                opttest = false;
                        }
                        double objval = dccgsub.model.GetObjValue();

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
                        multicutslist.Add(candidatesolset[ind]);

                        candidatesolset.RemoveAt(ind); costset.RemoveAt(ind);
                        counter++;
                    }
                }
                
                findw0rstscenario.Stop(); findw0rstscenariorecord = findw0rstscenario.ElapsedMilliseconds;
            }

            
            subproblemtime.Start();

            if (issubprofeasi)
            {
                
                
                if(omega < worst_scenario_cost)
                {
                    
                    IRange bendersoptimalitycut = ccgsub.GenBendersCut(Bendmaster, worst_u_j);

                    if (bendersoptimalitycut != null)
                    {
                        Add(bendersoptimalitycut);
                        Bendmaster.cutsStore.Add(bendersoptimalitycut);
                    }
                    subproblemtime.Stop();
                    subproblemtimerecord += subproblemtime.ElapsedMilliseconds;

                    multicutstime.Start();

                    
                    if (data.multiCutStrategy)
                    {
                        while (multicutslist.Count != 0)
                        {
                            int[] scenario = multicutslist[0];
                            multicutslist.RemoveAt(0);

                            ccgsub.ResetOBJaRHS_XI(y_solution, scenario);
                            bool solveflag = ccgsub.model.Solve();
                            if (solveflag)
                            {
                                IRange validcut = ccgsub.GenBendersCut(Bendmaster, scenario);
                                double subobj = ccgsub.model.GetObjValue();

                                if (validcut != null)
                                {
                                    Add(validcut);
                                    Bendmaster.cutsStore.Add(validcut);
                                }
                            }                            
                        }
                    }

                    multicutstime.Stop();
                    multicutstimerecord += multicutstime.ElapsedMilliseconds;

                    paretotime.Start();

                    
                    if (data.paretoCutStrategy)
                    {
                        BendersDualSub pareto = new BendersDualSub();
                        pareto.GenBendersubDual(y_solution, worst_u_j);
                        IRange paretoCut = pareto.GenParetoCut(y_solution, worst_u_j, worst_scenario_cost, data.pareto_y_sol, data.pareto_u_sol, Bendmaster);

                        if (paretoCut != null)
                        {
                            Add(paretoCut);
                            Bendmaster.cutsStore.Add(paretoCut);
                        }
                        pareto.model.End(); pareto = null;
                    }

                    paretotime.Stop();
                    paretotimerecord += paretotime.ElapsedMilliseconds;

                    if (updatePareto)
                    {
                        updatePareto = false;
                        data.pareto_y_sol = y_solution;
                        data.pareto_u_sol = worst_u_j;
                    }
                    
                }
                
            }
            else
            {
               
                IRange bendersfeasibilitycut = ccgsub.GenFeasibilityBendersCut(Bendmaster, y_solution);
                if (bendersfeasibilitycut != null)
                {
                    Add(bendersfeasibilitycut);
                    Bendmaster.cutsStore.Add(bendersfeasibilitycut);
                }

                subproblemtime.Stop();
                subproblemtimerecord += subproblemtime.ElapsedMilliseconds;
            }
            

            double bestbound = GetBestObjValue();
            Program.SolutionIteration.WriteLine($"{findw0rstscenariorecord / 1000}, {subproblemtimerecord / 1000.0}," +
                    $"{multicutstimerecord / 1000.0},{paretotimerecord / 1000.0},{worst_scenario_cost},{objectivecost},{bestbound},{omega}");
            
            Bendmaster.TreeDepth = Math.Max(Bendmaster.TreeDepth, GetCurrentNodeDepth());
            Bendmaster.NumofIter++;
            #endregion
        }        
    } 
    
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
        internal bool rootnodefathom;

        public SCNRLazyConsCallback(SCNRMaster Scenariomaster, Data data, int[] y_solution)
        {
            this.Scenariomaster = Scenariomaster;
            this.data = data;
            this.y_solution = y_solution;
            
            bestsubcost = -1;
            bestmastercost = -1;
            modifiedSub = new SCNRSub(y_solution, new int[data.DCSize], data);
            modifiedSub.GenScenarioSubproblem();
            solutionpool = new List<int[]>();
            searchtimes = new List<int>();
            rootnodefathom = true;
        }
        public override void Main()
        {
            #region
            int[] u_sol = new int[data.DCSize];
            double etaval = GetValue(Scenariomaster.eta);
            int numberofdisruptions = 0;
            for (int j = 0; j < data.DCSize; j++)
            {
                double val = GetValue(Scenariomaster.u_j[j]);
                if (val > 0.5)
                {
                    u_sol[j] = 1;
                    numberofdisruptions++;
                }
            }
            
            
            modifiedSub.ResetOBJaRHS_XI(y_solution, u_sol, data.big_M_r);
            bool issolved = modifiedSub.model.Solve();
            
            Cplex.Status solvingstatus = modifiedSub.model.GetStatus();

            if (issolved)
            {
                double benders_primal_cost = modifiedSub.model.GetObjValue();
                double[] x_r_values = new double[data.pathSize];                
                double[] z_values = modifiedSub.model.GetValues(modifiedSub.z_j);

                bool isfeaisble = true;
                for (int r = 0; r < data.pathSize; r++)
                {
                    x_r_values[r] = modifiedSub.model.GetValue(modifiedSub.x_r[r]);
                    
                    if (u_sol[data.s_r[r]] * x_r_values[r] >= data.epsilon)
                    {
                        if(x_r_values[r] < 0.001)
                        {
                            Program.Solvinglog.WriteLine("Error during finding the scenario using lazyconstraint callback: " + x_r_values[r]);
                        }
                                      
                        Scenariomaster.feasiblestatus = false;
                        Scenariomaster.bestFeasibleSolution = u_sol;

                        isfeaisble = false;

                        break;
                    }
                }

                if (isfeaisble)
                {                    
                    if (benders_primal_cost > Scenariomaster.bestObjVal)
                    {
                        Scenariomaster.bestObjVal = benders_primal_cost;
                        Scenariomaster.bestFeasibleSolution = u_sol;
                    }
                    else
                    {
                        bool sametobestsolution = true;
                        for (int j = 0; j < data.DCSize; j++)
                        {
                            if (Scenariomaster.bestFeasibleSolution[j] != u_sol[j])
                            {
                                sametobestsolution = false; break;
                            }
                        }

                    }

                    ILinearNumExpr extrem_point_exp = Scenariomaster.model.LinearNumExpr();
                                        
                    double subobj = 0;
                    for (int j = 0; j < data.DCSize; j++)
                    {
                        subobj += data.c_j[j]* z_values[j];

                    }
                    for (int r = 0; r < data.pathSize; r++)
                    {
                        subobj += data.d_r[r] * data.h_i[data.e_r[r]] * x_r_values[r];
                    }
                    for (int r = 0; r < data.pathSize; r++)
                    {
                        extrem_point_exp.AddTerm(- data.big_M_r[r] * x_r_values[r], Scenariomaster.u_j[data.s_r[r]]);
                    }

                    
                    extrem_point_exp.AddTerm(1, Scenariomaster.eta);
                    
                    Add(Scenariomaster.model.Le(extrem_point_exp, subobj));
                }                
            }
            else
            {
                Scenariomaster.bestFeasibleSolution = u_sol.ToArray();
                Scenariomaster.feasiblestatus = false;
                
            }

            if (!Scenariomaster.feasiblestatus)
                Abort();
            #endregion
        }
    }
    class SCNRMaster
    {
        #region
        public Cplex model;

        public INumVar eta;
        public INumVar[] u_j;
        Data data = new Data();
        public int number_of_var;
        public int number_of_con;

        public List<IRange> cutsStore;
        public List<IRange> optimalityCuts;
        public List<DataStructure> XiCoeff_list;
        public List<DataStructure> Temp_XiCoeff_list;
        
        
        public bool feasiblestatus = true;
        public int[] bestFeasibleSolution;
        public double bestObjVal = -float.MaxValue;
        public List<double> objvalpool = new List<double>();
        public List<int[]> solutionpool = new List<int[]>();

        internal int[] y_sol;
        public SCNRMaster(int[] y_sol)
        {
            this.y_sol = y_sol;
        }
        
        
        
        
        public void GenScenarioMaster()
        {
            #region
            model = new Cplex();

            u_j = new INumVar[data.DCSize];

            optimalityCuts = new List<IRange>();
            cutsStore = new List<IRange>();

            eta = model.NumVar(-float.MaxValue, float.MaxValue, NumVarType.Float);

            for (int j = 0; j < data.DCSize; j++)
            {
                
                u_j[j] = model.BoolVar($"u_{j}");
                number_of_var++;
            }
            model.AddMaximize(eta);

            int cons_1 = 1;
            if (cons_1 == 1)
            {
                ILinearNumExpr constraint = model.LinearNumExpr();               
                for (int j = 0; j < data.DCSize; j++)
                {
                    constraint.AddTerm(1, u_j[j]);
                }

                model.AddEq(constraint, data.max_dstroyed_DCs, "Disruption budget");                
                number_of_con++;
            }

            #endregion
        }
        public void Addpredefinedcuts(List<DataStructure> cutlist)
        {
            #region
            for (int c = 0; c < cutlist.Count; c++)
            {
                DataStructure dt = cutlist[c];
                ILinearNumExpr extrem_point_exp = model.LinearNumExpr();

                for (int r = 0; r < data.pathSize; r++)
                {
                    extrem_point_exp.AddTerm(-data.big_M_r[r] * dt.x_r_val[r], u_j[data.s_r[r]]);
                }
                
                extrem_point_exp.AddTerm(1, eta);
                
                model.AddLe(extrem_point_exp, dt.righthandsideU);
            }

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

        internal int[] u_sol;
        internal int[] y_sol;
        
        
        
        public SCNRSub(int[] y_sol, int[] u_sol, Data data)
        {
            this.y_sol = y_sol;            
            this.u_sol = u_sol;
            this.data = data;
        }
        
        
        
        
        public void GenScenarioSubproblem()
        {
            #region ccg subproblem model           
            

            model = new Cplex();            
                        
            
            
            number_of_var = 0; number_of_con = 0;
            
            x_r = new INumVar[data.pathSize];
            z_j = new INumVar[data.DCSize];            
            number_of_var = 0; number_of_con = 0;

            relatedtodual_V_value = new IRange[data.nodeSize];
            relatedtodual_W_value = new IRange[data.DCSize];
            relatedtodual_Delta_value = new IRange[data.pathSize];
            relatedtodual_Gamma_value = new IRange[data.pathSize];

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
                singleterm.AddTerm(data.d_r[r] * data.h_i[data.e_r[r]], x_r[r]);
            }
            for (int r = 0; r < data.pathSize; r++)
            {
                singleterm.AddTerm(data.big_M_r[r] * u_sol[data.s_r[r]], x_r[r]);
            }
            
            
            model.AddMinimize(singleterm);

            
            
            
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
                for (int l = 0; l < data.linkSize; l++)
                {
                    ILinearNumExpr consexpr = model.LinearNumExpr();
                    for (int r = 0; r < data.pathSize; r++)
                    {
                        if (data.pathlist[r].Contains(l))
                        {
                            consexpr.AddTerm(data.h_i[data.e_r[r]], x_r[r]);
                        }
                    }
                    model.AddLe(consexpr, data.linklist[l].LinkCapacity, $"section capacity_{l}");
                }
            }
            
            #endregion
        }
        
        
        
        
        
        public void ResetOBJaRHS_XI(int[] upd_y_solution, int[] upd_u_sol, double[] bigm_R)
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
                singleterm.AddTerm(data.d_r[r] * data.h_i[data.e_r[r]], x_r[r]);
            }
            for (int r = 0; r < data.pathSize; r++)
            {
                singleterm.AddTerm(bigm_R[r] * upd_u_sol[data.s_r[r]], x_r[r]);
            }
            
            
            mdsubobj.ClearExpr();
            
            mdsubobj.Expr = singleterm;

            
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
    class BendersMaster
    {
        public Cplex model;

        public INumVar omega;
        public INumVar[] y_j;
        
        public int number_of_var;
        public int number_of_con;
        public List<IRange> optimalityCuts;
        public List<IRange> feasibilityCuts;
        public List<IRange> cutsStore;

        public INumVar[] r_j;
        public INumVar[] u_j;

        public List<IConstraint> CCGoptimalityCuts;
        public int CCGfeasibilityCuts;
        public double L_bound = 0;

        List<IRange> generateCuts;
        public List<DataStructure> YCoeff_list;
        public List<DataStructure> Temp_YCoeff_list;

        Data data = new Data();
                
        public int NumofIter;
        public long TreeDepth;

        public void GenBDMasterproblem()
        {
            #region benders master model 

            model = new Cplex();

            number_of_var = 0; number_of_con = 0;
            y_j = new INumVar[data.DCSize];
            
            optimalityCuts = new List<IRange>();
            feasibilityCuts = new List<IRange>();
            cutsStore = new List<IRange>();

            CCGoptimalityCuts = new List<IConstraint>();

            CCGfeasibilityCuts = 0;

            omega = model.NumVar(0, float.MaxValue, NumVarType.Float, "omega");
            for (int j = 0; j < data.DCSize; j++)
            {
                y_j[j] = model.NumVar(0, 1, NumVarType.Int, $"y_{j}");
                number_of_var++; 
            }

            
            
            ILinearNumExpr numExpr = model.LinearNumExpr();

            for (int j = 0; j < data.DCSize; j++)
            {
                numExpr.AddTerm(data.f_j[j], y_j[j]);                
            }
            numExpr.AddTerm(1, omega);
            model.AddMinimize(numExpr);

            int cons_1 = 1;
            if (cons_1 == 1)
            {
                ILinearNumExpr totalOpenDCs = model.LinearNumExpr();
                for (int j = 0; j < data.DCSize; j++)
                {
                    totalOpenDCs.AddTerm(1, y_j[j]);
                }
                model.AddGe(totalOpenDCs, data.max_dstroyed_DCs + 1, $"the valid cut");
            }            
            #endregion
        }
        
        public void AddNewValidCuts(List<DataStructure> coeffList)
        {
            #region
            int[] u_sol = new int[data.DCSize];

            for (int d = 0; d < coeffList.Count; d++)
            {
                DataStructure dts = coeffList[d];
                ILinearNumExpr extrem_point_exp = model.LinearNumExpr();

                extrem_point_exp.AddTerm(1, omega);

                
                for (int r = 0; r < data.pathSize; r++)
                {                    
                    extrem_point_exp.AddTerm(- dts.delta_r_val[r] * (1 - u_sol[data.s_r[r]]), y_j[data.s_r[r]]);
                }

                
                double sumofvi = 0;
                for (int i = 0; i < data.nodeSize; i++)
                {
                    sumofvi += dts.v_i_val[i];
                }
                for (int a = 0; a < data.linkSize; a++)
                {                   
                    sumofvi += dts.pi_a_val[a] * data.linklist[a].LinkCapacity;
                }

                
                model.AddGe(extrem_point_exp, sumofvi);
            }
            #endregion
        }        
        
        public void GenLinearBDMasterproblem()
        {
            #region benders master model 

            model = new Cplex();
            YCoeff_list = new List<DataStructure>();
            Temp_YCoeff_list = new List<DataStructure>();
            generateCuts = new List<IRange>();

            number_of_var = 0; number_of_con = 0;
            y_j = new INumVar[data.DCSize];

            optimalityCuts = new List<IRange>();
            feasibilityCuts = new List<IRange>();
            cutsStore = new List<IRange>();

            CCGoptimalityCuts = new List<IConstraint>();

            CCGfeasibilityCuts = 0;

            omega = model.NumVar(0, float.MaxValue, NumVarType.Float, "omega");
            for (int j = 0; j < data.DCSize; j++)
            {
                y_j[j] = model.NumVar(0, 1, NumVarType.Float, $"y_{j}");
                number_of_var++;
            }

            
            
            INumExpr numExpr = model.NumExpr();

            for (int j = 0; j < data.DCSize; j++)
            {
                numExpr = model.Sum(model.Prod(data.f_j[j], y_j[j]), numExpr);
            }
            numExpr = model.Sum(omega, numExpr);
            model.AddMinimize(numExpr);

            int cons_1 = 1;
            if (cons_1 == 1)
            {
                ILinearNumExpr totalOpenDCs = model.LinearNumExpr();
                for (int j = 0; j < data.DCSize; j++)
                {
                    totalOpenDCs.AddTerm(1, y_j[j]);
                }
                model.AddGe(totalOpenDCs, data.max_dstroyed_DCs + 1, $"the valid cut");
            }
            int cons_2 = 2;
            if (cons_2 == 1)
            {
                ILinearNumExpr consExpr = model.LinearNumExpr();
                for (int j = 0; j < data.DCSize; j++)
                {
                    int corCommunityNode = Program.g_DCent_list[j].CorNodeSeqNo;
                    for (int l = 0; l < Program.g_node_list[corCommunityNode].OutgoingLinkList.Count; l++)
                    {
                        double arcdemand = Program.g_node_list[corCommunityNode].OutgoingLinkList[l].LinkCapacity;
                        consExpr.AddTerm(arcdemand, y_j[j]);
                    }
                }
                double alldemand = 0;
                for (int i = 0; i < data.nodeSize; i++)
                {
                    alldemand += data.h_i[i];
                }
                alldemand += (data.max_dstroyed_DCs + 1) * data.MinCap;

                model.AddGe(consExpr, alldemand);
            }
            #endregion
        }
        
        public void StabilizedBendersDec(double lambda, double eps_threshod, double alpha)
        {
            #region 
            double[] stab_y_sol = new double[data.DCSize];
            

            for (int j = 0; j < data.DCSize; j++)
            {
                stab_y_sol[j] = 1;
            }
            double mastercost = 0;
            double LB = -float.MaxValue;
            double benders_primal_cost = -float.MaxValue;
            double best_LP_bound = -float.MaxValue;
            double[] y_sol = new double[data.DCSize];

            int[] u_sol = new int[data.DCSize];
            bool modified = false;
            int iteration = 0;

            
            
            
            

            int It_counter = 0;
            
            
            

            
            bool masterfeasible = model.Solve();
            mastercost = model.GetObjValue();
            LB = mastercost;

            CCGSub ccgsub = new CCGSub();

            benders_primal_cost = -float.MaxValue;

            while (true)
            {
                double[] ast_y_sol = new double[data.DCSize];
                double[] one_matr = new double[data.DCSize];

                for (int j = 0; j < data.DCSize; j++)
                {
                    
                    ast_y_sol[j] = model.GetValue(y_j[j]);
                    one_matr[j] = 1;
                }

                
                

                if (It_counter >= data.Para_STB_IterLimit)
                {
                    
                    It_counter = 0;

                    
                    List<IRange> tempcutlist = new List<IRange>();
                    List<DataStructure> tempcoeff = new List<DataStructure>();
                    for (int c = 0; c < generateCuts.Count; c++)
                    {
                        double slackvalue = model.GetSlack(generateCuts[c]);
                        if (slackvalue < 0)
                        {
                            tempcutlist.Add(generateCuts[c]);
                            tempcoeff.Add(Temp_YCoeff_list[c]);
                        }
                        
                    }
                    generateCuts = tempcutlist;
                    Temp_YCoeff_list = tempcoeff;

                    if (lambda < 1)
                    {
                        lambda = 1;
                    }
                    else
                    {
                        if (lambda >= 1 && !modified)
                        {
                            eps_threshod = 0;
                            modified = true;
                        }
                        else
                        {
                            break;
                        }
                    }
                }

                if (data.Para_STB_strategy)
                {
                    
                    for (int j = 0; j < data.DCSize; j++)
                    {
                        stab_y_sol[j] = alpha * stab_y_sol[j] + (1 - alpha) * ast_y_sol[j];

                        y_sol[j] = lambda * ast_y_sol[j] + (1 - lambda) * stab_y_sol[j];
                    }
                }
                else
                {
                    
                    for (int j = 0; j < data.DCSize; j++)
                    {
                        stab_y_sol[j] = (ast_y_sol[j] + stab_y_sol[j]) / 2;

                        y_sol[j] = lambda * ast_y_sol[j] + (1 - lambda) * stab_y_sol[j] ;
                    }
                }

                Cplex.Status primalstatus = null;

                ccgsub = new CCGSub();
                ccgsub.GenlinearSubproblem(y_sol, u_sol);
                
                ccgsub.model.Solve();

                benders_primal_cost = ccgsub.model.GetObjValue();
                primalstatus = ccgsub.model.GetStatus();

                if (primalstatus == Cplex.Status.Optimal)
                {
                    DataStructure dts = new DataStructure(data);

                    double[] v_i_values = new double[data.nodeSize];
                    double[] delta_r_values = new double[data.pathSize];
                    double[] pi_a_values = new double[data.linkSize];

                    double sumofv = 0;
                    for (int i = 0; i < data.nodeSize; i++)
                    {
                        v_i_values[i] = ccgsub.model.GetDual(ccgsub.relatedtodual_V_value[i]);
                        sumofv += v_i_values[i];
                    }

                    for (int r = 0; r < data.pathSize; r++)
                    {
                        delta_r_values[r] = ccgsub.model.GetDual(ccgsub.relatedtodual_Delta_value[r]);
                    }

                    for (int a = 0; a < data.linkSize; a++)
                    {
                        
                        sumofv += pi_a_values[a] * data.linklist[a].LinkCapacity;
                    }

                    dts.pi_a_val = pi_a_values;
                    dts.delta_r_val = delta_r_values;
                    dts.v_i_val = v_i_values;

                    Temp_YCoeff_list.Add(dts);

                    ILinearNumExpr extrem_point_exp = model.LinearNumExpr();

                    extrem_point_exp.AddTerm(1, omega);

                    

                    for (int r = 0; r < data.pathSize; r++)
                    {
                        extrem_point_exp.AddTerm(-delta_r_values[r] * (1 - u_sol[data.s_r[r]]), y_j[data.s_r[r]]);
                    }

                    generateCuts.Add(model.AddGe(extrem_point_exp, sumofv));
                }
                else
                {
                    ILinearNumExpr combinatorialCut = model.LinearNumExpr();

                    double cardi = 0;
                    for (int j = 0; j < data.DCSize; j++)
                    {
                        if (y_sol[j] == 1)
                        {
                            combinatorialCut.AddTerm(-1, y_j[j]);
                            cardi++;
                        }
                        else
                        {
                            combinatorialCut.AddTerm(1, y_j[j]);
                        }
                    }
                    IRange cut = model.AddGe(combinatorialCut, 1 - cardi);
                }

                masterfeasible = model.Solve();
                mastercost = model.GetObjValue();
                LB = mastercost;

                if (LB > best_LP_bound)
                {
                    best_LP_bound = LB;
                    It_counter = 0;
                }
                else
                {
                    It_counter++;
                }
                iteration++;
            }
            #endregion
        }
    }
    class CCGMaster
    {
        public Cplex model;

        public INumVar omega;
        public INumVar[] y_j;
        
        public List<INumVar[]> x_l_r;        
        public List<INumVar[]> z_l_j;
        public List<int[]> u_l_j;

        public int number_of_var;
        public int number_of_con;
        public List<IRange> cutsStore;
        public List<IRange> feasiblecutsStore;

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
            u_l_j = new List<int[]>();
            cutsStore = new List<IRange>();
            feasiblecutsStore = new List<IRange>();

            for (int j = 0; j < data.DCSize; j++)
            {
                y_j[j] = model.BoolVar($"y_{j}");
                number_of_var++;
            }

            omega = model.NumVar(0, float.MaxValue, NumVarType.Float, "omega");

            ILinearNumExpr numExpr = model.LinearNumExpr();
            for (int j = 0; j < data.DCSize; j++)
            {
                numExpr.AddTerm(data.f_j[j], y_j[j]);                
            }
            numExpr.AddTerm(1, omega);            
            model.AddMinimize(numExpr);
            int cons_1 = 1;
            if(cons_1 == 1)
            {
                ILinearNumExpr consExpr = model.LinearNumExpr();
                for (int j = 0; j < data.DCSize; j++)
                {
                    consExpr.AddTerm(1, y_j[j]);
                }
                model.AddGe(consExpr, data.max_dstroyed_DCs + 1);
            }

            #endregion
        }
        
        public void GenCCGMasterproblem(int cutting_indicator)
        {
            #region model
            
            int[] worstscenario = u_l_j.Last();

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
            
            if (cutting_indicator == 1)
            {
                ILinearNumExpr singleterm = model.LinearNumExpr();
                for (int j = 0; j < data.DCSize; j++)
                {
                    singleterm.AddTerm(data.c_j[j], temp_z_j[j]);

                }
                for (int r = 0; r < data.pathSize; r++)
                {
                    singleterm.AddTerm(data.d_r[r]*data.h_i[data.e_r[r]], temp_x_r[r]);
                }                
                
                model.AddGe(omega, singleterm, "optimality cut");
            }

            
            int cons_1 = 1;
            if (cons_1 == 1)
            {
                for (int r = 0; r < data.pathSize; r++)
                {                    
                    model.AddLe(temp_x_r[r], model.Prod(y_j[data.s_r[r]], (1 - worstscenario[data.s_r[r]])), $"Transport_capacity_{r}_{data.s_r[r]}");
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
                for (int l = 0; l < data.linkSize; l++)
                {
                    ILinearNumExpr consexpr = model.LinearNumExpr();
                    for (int r = 0; r < data.pathSize; r++)
                    {
                        if (data.pathlist[r].Contains(l))
                        {
                            consexpr.AddTerm(data.h_i[data.e_r[r]], temp_x_r[r]);
                        }
                    }
                    model.AddLe(consexpr, data.linklist[l].LinkCapacity, $"section capacity_{l}");
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
        internal INumVar[,] f;
        internal INumVar[] s;
        internal INumVar[] l_j;


        public IRange[] relatedtodual_V_value;
        public IRange[] relatedtodual_W_value;
        public IRange[] relatedtodual_Delta_value;
        public IRange[] relatedtodual_Pi_value;

        public int number_of_var;
        public int number_of_con;
        
        
        
        
        
        public void GenCCGSubproblem(int[] y_sol, int[] u_sol)
        {
            #region ccg subproblem model            
            model = new Cplex();
            
            
           
            number_of_var = 0; number_of_con = 0;

            x_r = new INumVar[data.pathSize];
            z_j = new INumVar[data.DCSize];
            
            number_of_var = 0; number_of_con = 0;

            relatedtodual_V_value = new IRange[data.nodeSize];
            relatedtodual_W_value = new IRange[data.DCSize];
            relatedtodual_Delta_value = new IRange[data.pathSize];
            relatedtodual_Pi_value = new IRange[data.linkSize];

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
                singleterm.AddTerm(data.d_r[r] * data.h_i[data.e_r[r]], x_r[r]);
            }
            
            
            model.AddMinimize(singleterm);

            
            int cons_1 = 1;
            if (cons_1 == 1)
            {
                for (int r = 0; r < data.pathSize; r++)
                {
                    ILinearNumExpr cons = model.LinearNumExpr();
                    cons.AddTerm(1, x_r[r]);
                    IRange constraint = model.AddLe(cons,  y_sol[data.s_r[r]] * (1 - u_sol[data.s_r[r]]), $"Transport_capacity_{r}_{data.s_r[r]}");
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
                        if (data.e_r[r] == i)
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
                    contraint.AddTerm(1, z_j[j]);
                    
                    for (int r = 0; r < data.pathSize; r++)
                    {                        
                        if (data.s_r[r] == j)
                        {
                            contraint.AddTerm(-data.h_i[data.e_r[r]], x_r[r]);
                        }
                    }
                    
                    IRange constraint = model.AddGe(contraint, 0, $"purchasinggoods_{j}");
                    relatedtodual_W_value[j] = constraint;
                    number_of_con++;
                }
            }
            int cons_4 = 1;
            if(cons_4 == 1)
            {                
                for (int l = 0; l < data.linkSize; l++)
                {
                    ILinearNumExpr consexpr = model.LinearNumExpr();
                    for (int r = 0; r < data.pathSize; r++)
                    {
                        if (data.pathlist[r].Contains(l))
                        {
                            consexpr.AddTerm(data.h_i[data.e_r[r]], x_r[r]);
                        }                        
                    }
                    relatedtodual_Pi_value[l] = model.AddLe(consexpr, data.linklist[l].LinkCapacity, $"section capacity_{l}");
                }
            }
            #endregion
        }
        
        
        
        
        
        
        public void GenSimuCCGSub(int[] y_sol, int[] u_sol, double[] demand)
        {
            #region ccg subproblem model            
            model = new Cplex();
            
            

            number_of_var = 0; number_of_con = 0;

            x_r = new INumVar[data.pathSize];
            z_j = new INumVar[data.DCSize];

            number_of_var = 0; number_of_con = 0;

            relatedtodual_V_value = new IRange[data.nodeSize];
            relatedtodual_W_value = new IRange[data.DCSize];
            relatedtodual_Delta_value = new IRange[data.pathSize];
            relatedtodual_Pi_value = new IRange[data.linkSize];

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
                singleterm.AddTerm(data.d_r[r] * demand[data.e_r[r]], x_r[r]);
            }

            
            model.AddMinimize(singleterm);

            
            int cons_1 = 1;
            if (cons_1 == 1)
            {
                for (int r = 0; r < data.pathSize; r++)
                {
                    ILinearNumExpr cons = model.LinearNumExpr();
                    cons.AddTerm(1, x_r[r]);
                    IRange constraint = model.AddLe(cons, y_sol[data.s_r[r]] * (1 - u_sol[data.s_r[r]]), $"Transport_capacity_{r}_{data.s_r[r]}");
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
                        if (data.e_r[r] == i)
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
                    contraint.AddTerm(1, z_j[j]);

                    for (int r = 0; r < data.pathSize; r++)
                    {
                        if (data.s_r[r] == j)
                        {
                            contraint.AddTerm(-demand[data.e_r[r]], x_r[r]);
                        }
                    }

                    IRange constraint = model.AddGe(contraint, 0, $"purchasinggoods_{j}");
                    relatedtodual_W_value[j] = constraint;
                    number_of_con++;
                }
            }
            int cons_4 = 1;
            if (cons_4 == 1)
            {
                for (int l = 0; l < data.linkSize; l++)
                {
                    ILinearNumExpr consexpr = model.LinearNumExpr();
                    for (int r = 0; r < data.pathSize; r++)
                    {
                        if (data.pathlist[r].Contains(l))
                        {
                            consexpr.AddTerm(demand[data.e_r[r]], x_r[r]);
                        }

                    }
                    relatedtodual_Pi_value[l] = model.AddLe(consexpr, data.linklist[l].LinkCapacity, $"section capacity_{l}");
                }
            }
            #endregion
        }
        
        
        
        
        
        public void ResetOBJaRHS_XI(int[] upd_y_solution, int[] upd_u_sol)
        {
            #region

            
            for (int r = 0; r < data.pathSize; r++)
            {
                if (relatedtodual_Delta_value[r] != null)
                {
                    relatedtodual_Delta_value[r].UB = upd_y_solution[data.s_r[r]] * (1 - upd_u_sol[data.s_r[r]]);
                }
            }

            #endregion
        }
        
        
        
        
        
        
        public IRange GenBendersCut(BendersMaster BDM, int[] u_sol)
        {
            #region
            double[] v_i_values = new double[data.nodeSize];
            double[] delta_r_values = new double[data.pathSize];
            double[] pi_a_values = new double[data.linkSize];
            double[] w_j_values = new double[data.DCSize];

            double benders_primal_cost = model.GetObjValue();
            double sumofv = 0;
            for (int i = 0; i < data.nodeSize; i++)
            {
                v_i_values[i] = model.GetDual(relatedtodual_V_value[i]);
                sumofv += v_i_values[i];
            }
            for (int j = 0; j < data.DCSize; j++)
            {
                w_j_values[j] = model.GetDual(relatedtodual_W_value[j]);
            }
            for (int r = 0; r < data.pathSize; r++)
            {                
                delta_r_values[r] = model.GetDual(relatedtodual_Delta_value[r]);
            }

            for (int a = 0; a < data.linkSize; a++)
            {
                pi_a_values[a] = model.GetDual(relatedtodual_Pi_value[a]);                
                sumofv += pi_a_values[a] * data.linklist[a].LinkCapacity;
            }

            ILinearNumExpr extrem_point_exp = model.LinearNumExpr();

            extrem_point_exp.AddTerm(1, BDM.omega);

            
            for (int r = 0; r < data.pathSize; r++)
            {
                extrem_point_exp.AddTerm(- delta_r_values[r] * (1 - u_sol[data.s_r[r]]), BDM.y_j[data.s_r[r]]);
            }
            
            IRange cut = BDM.model.Ge(extrem_point_exp, sumofv);

            return cut;
            #endregion
        }
        
        
        
        
        
        
        
        public IRange GenFeasibilityBendersCut(BendersMaster BDM, int[] y_sol)
        {
            #region

            ILinearNumExpr combinatorialCut = BDM.model.LinearNumExpr();

            double cardi = 0;
            for (int j = 0; j < data.DCSize; j++)
            {
                if (y_sol[j] == 1)
                {
                    combinatorialCut.AddTerm(-1, BDM.y_j[j]);
                    cardi++;
                }
                else
                {
                    combinatorialCut.AddTerm(1, BDM.y_j[j]);
                }
            }
            IRange cut = BDM.model.Ge(combinatorialCut, 1 - cardi);
            return cut;

            #endregion
        }
        
        
        
        
        
        
        
        public void GenlinearSubproblem(double[] y_sol, int[] u_sol)
        {
            #region ccg subproblem model           

            model = new Cplex();
            
            

            number_of_var = 0; number_of_con = 0;

            x_r = new INumVar[data.pathSize];
            z_j = new INumVar[data.DCSize];

            number_of_var = 0; number_of_con = 0;

            relatedtodual_V_value = new IRange[data.nodeSize];
            relatedtodual_W_value = new IRange[data.DCSize];
            relatedtodual_Delta_value = new IRange[data.pathSize];
            relatedtodual_Pi_value = new IRange[data.linkSize];

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
                singleterm.AddTerm(data.d_r[r] * data.h_i[data.e_r[r]], x_r[r]);
            }

            
            model.AddMinimize(singleterm);

            
            int cons_1 = 1;
            if (cons_1 == 1)
            {
                for (int r = 0; r < data.pathSize; r++)
                {
                    ILinearNumExpr cons = model.LinearNumExpr();
                    cons.AddTerm(1, x_r[r]);
                    IRange constraint = model.AddLe(cons, y_sol[data.s_r[r]] * (1 - u_sol[data.s_r[r]]), $"Transport_capacity_{r}_{data.s_r[r]}");
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
                        if (data.e_r[r] == i)
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
                    contraint.AddTerm(1, z_j[j]);

                    for (int r = 0; r < data.pathSize; r++)
                    {
                        if (data.s_r[r] == j)
                        {
                            contraint.AddTerm(-data.h_i[data.e_r[r]], x_r[r]);
                        }
                    }

                    IRange constraint = model.AddGe(contraint, 0, $"purchasinggoods_{j}");
                    relatedtodual_W_value[j] = constraint;
                    number_of_con++;
                }
            }
            int cons_4 = 0;
            if (cons_4 == 1)
            {
                for (int l = 0; l < data.linkSize; l++)
                {
                    ILinearNumExpr consexpr = model.LinearNumExpr();
                    for (int r = 0; r < data.pathSize; r++)
                    {
                        if (data.pathlist[r].Contains(l))
                        {
                            consexpr.AddTerm(data.h_i[data.e_r[r]], x_r[r]);
                        }

                    }
                    relatedtodual_Pi_value[l] = model.AddLe(consexpr, data.linklist[l].LinkCapacity, $"section capacity_{l}");
                }
            }
            #endregion
        }
        
        public void ArcFLowModel(int[] y_sol, int[] u_sol)
        {
            #region
            model = new Cplex();
            
            

            number_of_var = 0; number_of_con = 0;
            f = new INumVar[data.DCSize, data.linkSize];
            s = new INumVar[data.DCSize];
            z_j = new INumVar[data.DCSize];
            l_j = new INumVar[data.DCSize];

            for (int j = 0; j < data.DCSize; j++)
            {
                z_j[j] = model.NumVar(0, double.MaxValue, NumVarType.Float, $"z_{j}");
                s[j] = model.NumVar(0, double.MaxValue, NumVarType.Float, $"s_{j}");
                l_j[j] = model.NumVar(0, double.MaxValue, NumVarType.Float, $"l_{j}");
                for (int a = 0; a < data.linkSize; a++)
                    f[j, a] = model.NumVar(0, double.MaxValue, NumVarType.Float, $"f_{j}_{a}");
            }

            
            ILinearNumExpr obj = model.LinearNumExpr();

            for (int j = 0; j < data.DCSize; j++)
            {
                obj.AddTerm(data.c_j[j], z_j[j]);
                for (int a = 0; a < data.linkSize; a++)
                    obj.AddTerm(data.linklist[a].LinkCost, f[j, a]);
            }

            model.AddMinimize(obj);

            
            
            for (int j = 0; j < data.DCSize; j++)
            {

                for (int i = 0; i < data.nodeSize; i++)
                {
                    ILinearNumExpr expr = model.LinearNumExpr();
                    foreach (Link a in Program.g_node_list[i].OutgoingLinkList)
                        expr.AddTerm(1.0, f[j, a.LinkSeqNo]);

                    model.AddLe(expr, data.big_M * y_sol[j] * (1 - u_sol[j]), $"activation_{j}");
                }

            }
            
            for (int j = 0; j < data.DCSize; j++)
            {
                int nodeseq = Program.g_DCent_list[j].CorNodeSeqNo;

                ILinearNumExpr expr = model.LinearNumExpr();

                foreach (Link a in Program.g_node_list[nodeseq].OutgoingLinkList)
                    expr.AddTerm(1.0, f[j, a.LinkSeqNo]);

                foreach (Link a in Program.g_node_list[nodeseq].IngoingLinkList)
                    expr.AddTerm(-1.0, f[j, a.LinkSeqNo]);

                model.AddEq(s[j], expr, $"commodity-specific_{j}");
            }
            

            for (int j = 0; j < data.DCSize; j++)
            {

                ILinearNumExpr rightexpr = model.LinearNumExpr();
                rightexpr.AddTerm(y_sol[j] * (1 - u_sol[j]), z_j[j]);

                model.AddLe(model.Sum(s[j], l_j[j]), rightexpr, $"SupplyInjection_{j}");

                
            }

            for (int i = 0; i < data.nodeSize; i++)
            {
                ILinearNumExpr expr = model.LinearNumExpr();

                for (int j = 0; j < data.DCSize; j++)
                {
                    if (Program.g_node_list[i].CorDcenterSeqNo == j)
                        continue;
                    
                    foreach (Link a in Program.g_node_list[i].IngoingLinkList)
                        expr.AddTerm(1.0, f[j, a.LinkSeqNo]);

                    foreach (Link a in Program.g_node_list[i].OutgoingLinkList)
                        expr.AddTerm(-1.0, f[j, a.LinkSeqNo]);

                }

                
                if (Program.g_node_list[i].DCstate != 0)
                {
                    int dcno = Program.g_node_list[i].CorDcenterSeqNo;
                    
                    if (y_sol[dcno] * (1 - u_sol[dcno]) == 1)
                        expr.AddTerm(1.0, l_j[dcno]);
                }

                model.AddGe(expr, data.h_i[i], $"Demand_{i}");
            }

            for (int a = 0; a < data.linkSize; a++)
            {
                ILinearNumExpr expr = model.LinearNumExpr();
                for (int j = 0; j < data.DCSize; j++)
                    expr.AddTerm(1.0, f[j, a]);

                model.AddLe(expr, data.linklist[a].LinkCapacity, $"ArcCap_{a}");
            }
            model.ExportModel("NFMmodel.lp");
            Console.WriteLine();
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
        public INumVar[] u_j;
        public INumVar[] delta_r;
        public INumVar[] b_r;
        public INumVar[] pi_a;
        public IRange[] relatedtoUB_delta_value;
        public IRange[] relatedto_DualX_value;

        public int number_of_var;
        public int number_of_con;
        
        public void GenDualofCCGsub(int[] y_sol, double[] UB_delta)
        {
            #region
            model = new Cplex();
            number_of_var = 0; number_of_con = 0;
            
            v_i = new INumVar[data.nodeSize];
            w_j = new INumVar[data.DCSize];
            u_j = new INumVar[data.DCSize];
            pi_a = new INumVar[data.linkSize];
            delta_r = new INumVar[data.pathSize];
            b_r = new INumVar[data.pathSize];
            relatedtoUB_delta_value = new IRange[data.pathSize];
            relatedto_DualX_value = new IRange[data.pathSize];

            for (int r = 0; r < data.pathSize; r++)
            {
                delta_r[r] = model.NumVar(float.MinValue, 0, NumVarType.Float, $"delta_{r}");
                number_of_var++;
                b_r[r] = model.NumVar(float.MinValue, 0, NumVarType.Float, $"b_{r}");
                number_of_var++;
            }
            for (int i = 0; i < data.nodeSize; i++)
            {
                v_i[i] = model.NumVar(0, float.MaxValue, NumVarType.Float, $"v_{i}");
                number_of_var++;
            }
            for (int j = 0; j < data.DCSize; j++)
            {
                w_j[j] = model.NumVar(float.MinValue, float.MaxValue, NumVarType.Float, $"w_{j}");
                number_of_var++;
                u_j[j] = model.NumVar(0, y_sol[j], NumVarType.Int, $"u_{j}");
                number_of_var++;
            }
            for (int a = 0; a < data.linkSize; a++)
            {
                pi_a[a] = model.NumVar(float.MinValue, 0, NumVarType.Float, $"pi_{a}");
                number_of_var++;
            }

            
            ILinearNumExpr subobj = model.LinearNumExpr();
            for (int r = 0; r < data.pathSize; r++)
            {
                subobj.AddTerm(y_sol[data.s_r[r]], delta_r[r]);
            }
            for (int r = 0; r < data.pathSize; r++)
            {
                subobj.AddTerm(-y_sol[data.s_r[r]], b_r[r]);
            }
            for (int i = 0; i < data.nodeSize; i++)
            {
                subobj.AddTerm(1, v_i[i]);
            }            

            for (int a = 0; a < data.linkSize; a++)
            {
                subobj.AddTerm(data.linklist[a].LinkCapacity, pi_a[a]);
            }

            model.AddMaximize(subobj);

            int cons_0 = 1;
            if (cons_0 == 1)
            {                
                for (int r = 0; r < data.pathSize; r++)
                {
                    ILinearNumExpr consexpr = model.LinearNumExpr();
                    consexpr.AddTerm(data.h_i[data.e_r[r]], w_j[data.s_r[r]]);

                    consexpr.AddTerm(1, v_i[data.e_r[r]]);

                    consexpr.AddTerm(1, delta_r[r]);

                    for (int l = 0; l < data.linkSize; l++)
                    {
                        if (data.pathlist[r].Contains(l))
                        {
                            consexpr.AddTerm(data.h_i[data.e_r[r]], pi_a[l]);
                        }
                    }

                    model.AddLe(consexpr, data.d_r[r] * data.h_i[data.e_r[r]], $"X_dual[{r}]");
                }
            }

            int cons_1 = 1;
            if (cons_1 == 1)
            {
                for (int j = 0; j < data.DCSize; j++)
                {
                    IRange constraint = model.AddEq(w_j[j], -data.c_j[j], $"Z_Dual_{j}");
                    number_of_con++;
                }
            }
            int cons_2 = 1;
            if (cons_2 == 1)
            {                
                for (int r = 0; r < data.pathSize; r++)
                {
                    ILinearNumExpr consexpr = model.LinearNumExpr();
                    consexpr.AddTerm(1, b_r[r]);
                    consexpr.AddTerm(-UB_delta[r], u_j[data.s_r[r]]);

                    IRange constraint = model.AddGe(consexpr, 0, $"BigM_dual_{r}");
                    relatedtoUB_delta_value[r] = constraint;
                }                
            }
            int cons_3 = 1;
            if (cons_3 == 1)
            {
                for (int r = 0; r < data.pathSize; r++)
                {
                    model.AddGe(b_r[r], delta_r[r], $"BigM2_dual_{r}");
                    
                }
            }            
            int cons_4 = 1;
            if(cons_4 == 1)
            {
                ILinearNumExpr consexpr = model.LinearNumExpr();
                for (int j = 0; j < data.DCSize; j++)
                {
                    consexpr.AddTerm(1, u_j[j]);
                }
                model.AddEq(consexpr, data.max_dstroyed_DCs, "budget uncertainty");
            }
            
            #endregion
        }
        
        public void ResetdualccgSubObjcons(int[] upd_y_sol, double[] upd_UB_delta)
        {
            #region
            
            for (int j = 0; j < data.DCSize; j++)
            {
                u_j[j].UB = upd_y_sol[j];
            }
            
            IObjective originalobj = model.GetObjective();
            originalobj.ClearExpr();

            
            ILinearNumExpr subobj = model.LinearNumExpr();
            for (int r = 0; r < data.pathSize; r++)
            {
                subobj.AddTerm(upd_y_sol[data.s_r[r]], delta_r[r]);
            }
            for (int r = 0; r < data.pathSize; r++)
            {
                subobj.AddTerm(-upd_y_sol[data.s_r[r]], b_r[r]);
            }
            for (int i = 0; i < data.nodeSize; i++)
            {
                subobj.AddTerm(1, v_i[i]);
            }
            for (int a = 0; a < data.linkSize; a++)
            {
                subobj.AddTerm(data.linklist[a].LinkCapacity, pi_a[a]);
            }
            originalobj.Expr = subobj;

            
            
            for (int r = 0; r < data.pathSize; r++)
            {
                
                relatedtoUB_delta_value[r].ClearExpr();
                
                ILinearNumExpr consexpr = model.LinearNumExpr();
                consexpr.AddTerm(1, b_r[r]);
                consexpr.AddTerm(-upd_UB_delta[r], u_j[data.s_r[r]]);
                
                relatedtoUB_delta_value[r].Expr = consexpr;
            }
            #endregion
        }
        
    }
    class BendersDualSub
    {
        public Cplex model;
        public Data data = new Data();

        public INumVar[] v_i;
        public INumVar[] delta_r;        
        public INumVar[] w_j;
        public INumVar[] pi_a;
        
        public int number_of_var;
        public int number_of_con;
        public List<IConstraint> relatedto_X_dual;
        public IConstraint[,] relatedto_B_dual;

        public void GenBendersubDual(int[] y_sol, int[] u_sol)
        {
            #region
            model = new Cplex();
            number_of_var = 0; number_of_con = 0;
            
            v_i = new INumVar[data.nodeSize];
            w_j = new INumVar[data.DCSize];            
            delta_r = new INumVar[data.pathSize];
            pi_a = new INumVar[data.linkSize];
            INumVar[] gamma_r = new INumVar[data.pathSize];
            for (int r = 0; r < data.pathSize; r++)
            {
                delta_r[r] = model.NumVar(float.MinValue, 0, NumVarType.Float, $"delta_{r}");
                number_of_var++;               
            }
            for (int i = 0; i < data.nodeSize; i++)
            {
                v_i[i] = model.NumVar(0, float.MaxValue, NumVarType.Float, $"v_{i}");
                number_of_var++;
            }
            for (int j = 0; j < data.DCSize; j++)
            {
                w_j[j] = model.NumVar(float.MinValue, float.MaxValue, NumVarType.Float, $"w_{j}");
                number_of_var++;                
            }
            for (int a = 0; a < data.linkSize; a++)
            {
                pi_a[a] = model.NumVar(float.MinValue, 0, NumVarType.Float, $"pi_{a}");
                number_of_var++;
            }

            
            ILinearNumExpr subobj = model.LinearNumExpr();
            for (int r = 0; r < data.pathSize; r++)
            {
                subobj.AddTerm(y_sol[data.s_r[r]]*(1 - u_sol[data.s_r[r]]), delta_r[r]);
                
            }
            for (int i = 0; i < data.nodeSize; i++)
            {
                subobj.AddTerm(1, v_i[i]);
            }            

            for (int a = 0; a < data.linkSize; a++)
            {
                subobj.AddTerm(data.linklist[a].LinkCapacity, pi_a[a]);
            }
            model.AddMaximize(subobj);

            int cons_0 = 1;
            if (cons_0 == 1)
            {                
                for (int r = 0; r < data.pathSize; r++)
                {
                    ILinearNumExpr consexpr = model.LinearNumExpr();
                    consexpr.AddTerm(data.h_i[data.e_r[r]], w_j[data.s_r[r]]);

                    consexpr.AddTerm(1, v_i[data.e_r[r]]);

                    consexpr.AddTerm(1, delta_r[r]);

                    for (int l = 0; l < data.linkSize; l++)
                    {
                        if (data.pathlist[r].Contains(l))
                        {
                            consexpr.AddTerm(data.h_i[data.e_r[r]], pi_a[l]);
                        }
                    }

                    model.AddLe(consexpr, data.d_r[r] * data.h_i[data.e_r[r]], $"X_dual[{r}]");
                }
            }

            int cons_1 = 1;
            if (cons_1 == 1)
            {
                for (int j = 0; j < data.DCSize; j++)
                {
                    IRange constraint = model.AddEq(w_j[j], - data.c_j[j], $"Z_Dual_{j}");
                    number_of_con++;
                }
            }            
            int cons_2 = 0;
            if (cons_2 == 1)
            {
                for (int i = 0; i < data.nodeSize; i++)
                {
                    model.AddLe(v_i[i], data.B_i[i], $"upper bound of vi_{i}");
                }
            }
            #endregion
        }
        
        public void GendeltaUBsubDual(int[] y_sol, int[] u_sol)
        {
            #region
            model = new Cplex();
            number_of_var = 0; number_of_con = 0;
            
            v_i = new INumVar[data.nodeSize];
            w_j = new INumVar[data.DCSize];
            delta_r = new INumVar[data.pathSize];
            
            for (int r = 0; r < data.pathSize; r++)
            {
                delta_r[r] = model.NumVar(float.MinValue, 0, NumVarType.Float, $"delta_{r}");
                number_of_var++;
            }
            for (int i = 0; i < data.nodeSize; i++)
            {
                v_i[i] = model.NumVar(0, float.MaxValue, NumVarType.Float, $"v_{i}");
                number_of_var++;
            }
            for (int j = 0; j < data.DCSize; j++)
            {
                w_j[j] = model.NumVar(float.MinValue, float.MaxValue, NumVarType.Float, $"w_{j}");
                number_of_var++;
            }
            

            
            ILinearNumExpr subobj = model.LinearNumExpr();
            for (int r = 0; r < data.pathSize; r++)
            {
                subobj.AddTerm(y_sol[data.s_r[r]] * (1 - u_sol[data.s_r[r]]), delta_r[r]);
            }
            for (int i = 0; i < data.nodeSize; i++)
            {
                subobj.AddTerm(1, v_i[i]);
            }           
            
            model.AddMaximize(subobj);

            int cons_0 = 1;
            if (cons_0 == 1)
            {

                for (int r = 0; r < data.pathSize; r++)
                {
                    ILinearNumExpr consexpr = model.LinearNumExpr();
                    consexpr.AddTerm(data.h_i[data.e_r[r]], w_j[data.s_r[r]]);

                    consexpr.AddTerm(1, v_i[data.e_r[r]]);

                    consexpr.AddTerm(1, delta_r[r]);

                    model.AddLe(consexpr, data.d_r[r] * data.h_i[data.e_r[r]], $"X_dual[{r}]");
                }
            }

            int cons_1 = 1;
            if (cons_1 == 1)
            {
                for (int j = 0; j < data.DCSize; j++)
                {
                    IRange constraint = model.AddEq(w_j[j], -data.c_j[j], $"Z_Dual_{j}");
                    number_of_con++;
                }
            }
            int cons_2 = 1;
            if (cons_2 == 1)
            {
                for (int i = 0; i < data.nodeSize; i++)
                {
                    model.AddLe(v_i[i], data.B_i[i], $"upper bound of vi_{i}");
                }
            }
            #endregion
        }
        
        public IRange GenParetoCut(int[] y_sol, int[] u_sol, double objcostofdualsub,
            int[] y0_sol, int[] u0_sol, BendersMaster BDM)
        {
            #region
            IObjective original_objf = model.GetObjective();

            int objflag = 1;
            if (objflag == 1)
            {
                ILinearNumExpr subobj = model.LinearNumExpr();
                for (int r = 0; r < data.pathSize; r++)
                {
                    subobj.AddTerm(y0_sol[data.s_r[r]] * (1 - u0_sol[data.s_r[r]]), delta_r[r]);
                }
                for (int i = 0; i < data.nodeSize; i++)
                {
                    subobj.AddTerm(1, v_i[i]);
                }
                for (int a = 0; a < data.linkSize; a++)
                {
                    subobj.AddTerm(data.linklist[a].LinkCapacity, pi_a[a]);
                }
                original_objf.ClearExpr();
                original_objf.Expr = subobj;

            }
            
            int cons = 1;
            if (cons == 1)
            {
                ILinearNumExpr subobj = model.LinearNumExpr();
                for (int r = 0; r < data.pathSize; r++)
                {
                    subobj.AddTerm(y_sol[data.s_r[r]] * (1 - u_sol[data.s_r[r]]), delta_r[r]);
                }
                for (int i = 0; i < data.nodeSize; i++)
                {
                    subobj.AddTerm(1, v_i[i]);
                }
                for (int a = 0; a < data.linkSize; a++)
                {
                    subobj.AddTerm(data.linklist[a].LinkCapacity, pi_a[a]);
                }
                model.AddEq(subobj, objcostofdualsub, "Pareto optimal cut");
            }

            try
            {
                model.Solve();
                Console.WriteLine(model.GetStatus());
                double objcost = model.GetObjValue();

                if (objcost > objcostofdualsub)
                {
                    double[] v_i_values = new double[data.nodeSize];
                    double[] delta_r_values = new double[data.pathSize];
                    double[] pi_a_values = new double[data.linkSize];

                    double benders_primal_cost = model.GetObjValue();
                    double sumofv = 0;
                    for (int i = 0; i < data.nodeSize; i++)
                    {
                        v_i_values[i] = model.GetValue(v_i[i]);
                        sumofv += v_i_values[i];
                    }

                    for (int r = 0; r < data.pathSize; r++)
                    {
                        delta_r_values[r] = model.GetValue(delta_r[r]);
                    }
                    for (int a = 0; a < data.linkSize; a++)
                    {
                        pi_a_values[a] = model.GetValue(pi_a[a]);
                        sumofv += pi_a_values[a] * data.linklist[a].LinkCapacity;
                    }
                    ILinearNumExpr extrem_point_exp = model.LinearNumExpr();

                    extrem_point_exp.AddTerm(1, BDM.omega);

                    
                    for (int r = 0; r < data.pathSize; r++)
                    {
                        extrem_point_exp.AddTerm(-delta_r_values[r] * (1 - u_sol[data.s_r[r]]), BDM.y_j[data.s_r[r]]);
                    }

                    IRange cut = BDM.model.Ge(extrem_point_exp, sumofv);
                    return cut;
                }
            }
            catch (ILOG.Concert.Exception e)
            {
                Program.Solvinglog.WriteLine("Error during finding the pareto cut in solving subproblem: " + e.Message);

            }
            return null;
            #endregion
        }
        
        public void ResetBDDSubObjforDelta(int[] upd_y_sol, int[] upd_u_j)
        {
            #region
            IObjective ccgsubobj = model.GetObjective();
            ccgsubobj.ClearExpr();
            ILinearNumExpr subobj = model.LinearNumExpr();
            for (int r = 0; r < data.pathSize; r++)
            {
                subobj.AddTerm(upd_y_sol[data.s_r[r]] * (1 - upd_u_j[data.s_r[r]]), delta_r[r]);
            }
            for (int i = 0; i < data.nodeSize; i++)
            {
                subobj.AddTerm(1, v_i[i]);
            }            
            ccgsubobj.Expr = subobj;

            #endregion
        }
    }
    class CACGBD
    {
        static Data data = new Data();
        Solution solution = new Solution();

        public List<int[]> xI_sols_pool;
        public double scenarioCost;

        public SCNRMaster BranchandCutForWorstCsenario(int[] y_sol)
        {
            #region            
            SCNRMaster ssmaster = new SCNRMaster(y_sol);
            xI_sols_pool = new List<int[]>();
            Random rand = new Random();

            
            ssmaster.GenScenarioMaster();

            List<DataStructure> p_predefinedcut_list = new List<DataStructure>();

            ssmaster.Addpredefinedcuts(p_predefinedcut_list);

            
            ssmaster.model.SetParam(Cplex.Param.MIP.Strategy.Search, Cplex.MIPSearch.Traditional);
            ssmaster.model.SetParam(Cplex.Param.Preprocessing.Presolve, false);
            ssmaster.model.SetParam(Cplex.Param.Threads, 1);

            
            

            
            ssmaster.model.SetParam(Cplex.Param.MIP.Pool.Capacity, data.maxMultiCuts); 
            ssmaster.model.SetParam(Cplex.Param.MIP.Pool.Replace, 2); 
                                                                      
                                                                      
                                                                      

            ssmaster.model.SetParam(Cplex.Param.MIP.Pool.Intensity, 0); 
                                                                        
                                                                        
                                                                        

            ssmaster.model.Use(new SCNRLazyConsCallback (ssmaster, data, y_sol));

            ssmaster.model.Solve();

            if (!ssmaster.feasiblestatus)
                return ssmaster;

            Cplex.Status modelsts = ssmaster.model.GetStatus();
            scenarioCost = ssmaster.model.GetObjValue();
            double bestScenarioCost = ssmaster.model.GetBestObjValue();

            Console.WriteLine($"Optimal objective cost: {scenarioCost}");
            Cplex.Status status = ssmaster.model.GetStatus();
            int[] u_sol = new int[data.DCSize];
            for (int j = 0; j < data.DCSize; j++)
            {
                double val = ssmaster.model.GetValue(ssmaster.u_j[j]);
                if (val > 0.5)
                {
                    u_sol[j] = 1;
                }
            }

            
            if (ssmaster.bestObjVal - scenarioCost > 1e-4)
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
                    int[] feasiblesolution = new int[data.DCSize];
                    for (int j = 0; j < data.DCSize; j++)
                    {
                        if (ssmaster.model.GetValue(ssmaster.u_j[j], n) > 0.5)
                        {
                            feasiblesolution[j] = 1;
                        }
                        if (feasiblesolution[j] != u_sol[j])
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
        public void generateParetoPoint()
        {
            #region
            CCGMaster ccgmst = new CCGMaster();
            
            ccgmst.InitializeCCGMaster();
            int[] u_sol = new int[data.DCSize];
            
            
            
            
            
            

            
            
            

            data.pareto_y_sol = new int[data.DCSize];
            
            for (int j = 0; j < data.DCSize; j++)
            {
                
                
                
                
                
                data.pareto_y_sol[j] = 1;
            }
            ccgmst.model.End();
            ccgmst = null;

            SCNRMaster submaster = BranchandCutForWorstCsenario(data.pareto_y_sol);

            Cplex.Status cur_status = submaster.model.GetStatus();
            data.pareto_u_sol = new int[data.DCSize];
            for (int l = 0; l < data.DCSize; l++)
            {
                var u_val = submaster.model.GetValue(submaster.u_j[l]);

                if (u_val > 0.5)
                {
                    data.pareto_u_sol[l] = 1;
                }
            }
            submaster = null;
            #endregion
        }        
        
        public BendersMaster BranchandCut()
        {
            #region
            Program.SolutionIteration.WriteLine("findworsttime, subproblemtime,multicutstime,paretotime,worst_scenario_cost, objective cost, bestbound, Omegaval");

            Program.g_comparativeIndicators.WriteLine("# feasol, # optsol, best integer, best bound, Gap, CPU time");
            BendersMaster bdmaster = new BendersMaster();
            
            bdmaster.GenBDMasterproblem();

            BendersMaster linearbdm = new BendersMaster();

            CCGMaster CCGMasterProblem = new CCGMaster();
            
            
            List<DataStructure> YCoeff_list = new List<DataStructure>();
            
            Stopwatch optimize_BAC = new Stopwatch();
            CCGMasterProblem.YCoeff_list = new List<DataStructure>();

            if (data.stablilization_at_rootnode)
            {
                int count = 0;
                while (true)
                {                    
                    
                    int total_sets = 3;
                    double[] alpha_list = new double[total_sets];
                    double[] lambda_list = new double[total_sets];
                    
                    for (int m = 0; m < total_sets; m++)
                    {
                        double lambda_val = data.Para_STB_Lambada + 0.1 * m;
                        double alpha_val = data.Para_STB_Alpha - 0.1 * m;

                        alpha_list[m] = alpha_val;
                        lambda_list[m] = lambda_val;
                    }

                    object lockObject = new object();
                    var options = new ParallelOptions() { MaxDegreeOfParallelism = Environment.ProcessorCount };


                    Parallel.For(0, total_sets, options, (i, state) =>
                    {
                        BendersMaster bdmp = new BendersMaster();
                        bdmp.GenLinearBDMasterproblem();
                        bdmp.AddNewValidCuts(YCoeff_list);
                        bdmp.StabilizedBendersDec(lambda_list[i], data.Para_STB_Epsilon, alpha_list[i]);

                        lock (lockObject)
                        {
                            for (int f = 0; f < bdmp.Temp_YCoeff_list.Count; f++)
                            {
                                YCoeff_list.Add(bdmp.Temp_YCoeff_list[f]);
                            }
                        }
                    });

                    count++;
                    
                    if (count >= data.rootnodesolverounds)
                    {
                        break;
                    }                    
                }                
            }

            if (data.paretoCutStrategy)
            {
                generateParetoPoint();
            }

            
            bdmaster.AddNewValidCuts(YCoeff_list);
            CCGMasterProblem = null; 
            
            bdmaster.model.SetParam(Cplex.Param.TimeLimit, data.TL);
            
            bdmaster.model.SetParam(Cplex.Param.Preprocessing.Presolve, false);
            bdmaster.model.SetParam(Cplex.Param.Threads, 1);
            
            

            bdmaster.model.SetParam(Cplex.Param.MIP.Strategy.Search, Cplex.MIPSearch.Traditional);
            bdmaster.model.Use(new BendersLazyConsCallback(bdmaster, data));
            
            bdmaster.model.SetOut(Program.TWoutput);
            int[] y_solution = new int[data.DCSize];
            double omegaval = 0;
            try
            {
                optimize_BAC.Start();
                bdmaster.model.Solve();
                optimize_BAC.Stop();

                double optobj = bdmaster.model.GetObjValue();
                double bestbound = bdmaster.model.GetBestObjValue();
                double totalCPUtime = optimize_BAC.ElapsedMilliseconds / 1000;
                double relativeGap = bdmaster.model.GetMIPRelativeGap();
                int totalSolutionCount = bdmaster.model.GetSolnPoolNsolns();
                int Nofeasols = 0;
                int Nooptsols = 0;
                
                for (int j = 0; j < data.DCSize; j++)
                {
                    var yval = bdmaster.model.GetValue(bdmaster.y_j[j]);
                    if(yval > 0.5)
                    {
                        y_solution[j] = 1;
                    }                     
                }
                omegaval = bdmaster.model.GetValue(bdmaster.omega);

                List<double> optlist = new List<double>();
                List<double> feasilist = new List<double>();

                if (totalSolutionCount >= 1)
                {
                    for (int n = 0; n < totalSolutionCount; n++)
                    {
                        int[] y_val = new int[data.DCSize];
                        double omega_val = 0;

                        for (int j = 0; j < data.DCSize; j++)
                        {
                            if (bdmaster.model.GetValue(bdmaster.y_j[j], n) > 0.5)
                            {
                                y_val[j] = 1;
                            }
                            omega_val = bdmaster.model.GetValue(bdmaster.omega, n);
                        }
                        double currentobj = 0;
                        for (int j = 0; j < data.DCSize; j++)
                        {
                            currentobj += data.f_j[j] * y_val[j];
                        }
                        currentobj += omega_val;

                        if (currentobj > optobj)
                        {
                            Nofeasols++;
                            feasilist.Add(currentobj);
                        }
                        else
                        {
                            optlist.Add(currentobj);
                            Nooptsols++;
                        }
                    }
                }

                Program.g_comparativeIndicators.WriteLine($"{Nofeasols}, {Nooptsols}, {optobj}, {bestbound}, {relativeGap}, {totalCPUtime}");
                Program.g_comparativeIndicators.WriteLine();

                Program.g_comparativeIndicators.Write("feasiblesol" + ",");

                for (int i = 0; i < Nofeasols; i++)
                {
                    Program.g_comparativeIndicators.Write($"{feasilist[i]}" + ",");
                }
                Program.g_comparativeIndicators.WriteLine();

                Program.g_comparativeIndicators.Write("Optimalsol" + ",");

                for (int i = 0; i < Nooptsols; i++)
                {
                    Program.g_comparativeIndicators.Write($"{optlist[i]}" + ",");
                }

                Program.SolutionIteration.WriteLine(); 
                
                int NumofNodesExplored = bdmaster.model.Nnodes;

                Program.g_comparativeIndicators.WriteLine();
                Program.g_comparativeIndicators.WriteLine("TreeDepth, # NodesExplored, #Iter, Obj,bestbound, Gap, CPU time");
                Program.g_comparativeIndicators.WriteLine($"{bdmaster.TreeDepth},{NumofNodesExplored},{bdmaster.NumofIter}, {optobj},{bestbound}, {relativeGap}, {totalCPUtime}");

                Cplex.Status solvests = bdmaster.model.GetStatus();

                Program.SolutionIteration.WriteLine($"Total B&C time,{totalCPUtime}");
                Program.SolutionIteration.WriteLine($"Optimal obj,{optobj}");
                Program.SolutionIteration.WriteLine($"Best Linear bound,{bestbound}");

            }
            catch (ILOG.Concert.Exception e)
            {

                double optobj = bdmaster.model.GetObjValue();
                double bestbound = bdmaster.model.GetBestObjValue();
                double totalCPUtime = optimize_BAC.ElapsedMilliseconds / 1000;
                double relativeGap = bdmaster.model.GetMIPRelativeGap();
                int totalSolutionCount = bdmaster.model.GetSolnPoolNsolns();
                int Nofeasols = 0;
                int Nooptsols = 0;

                for (int j = 0; j < data.DCSize; j++)
                {
                    var yval = bdmaster.model.GetValue(bdmaster.y_j[j]);
                    if (yval > 0.5)
                    {
                        y_solution[j] = 1;
                    }
                }
                omegaval = bdmaster.model.GetValue(bdmaster.omega);

                List<double> optlist = new List<double>();
                List<double> feasilist = new List<double>();

                if (totalSolutionCount >= 1)
                {
                    for (int n = 0; n < totalSolutionCount; n++)
                    {
                        int[] y_val = new int[data.DCSize];
                        double omega_val = 0;

                        for (int j = 0; j < data.DCSize; j++)
                        {
                            if (bdmaster.model.GetValue(bdmaster.y_j[j], n) > 0.5)
                            {
                                y_val[j] = 1;
                            }
                            omega_val = bdmaster.model.GetValue(bdmaster.omega, n);
                        }
                        double currentobj = 0;
                        for (int j = 0; j < data.DCSize; j++)
                        {
                            currentobj += data.f_j[j] * y_val[j];
                        }
                        currentobj += omega_val;

                        if (currentobj > optobj)
                        {
                            Nofeasols++;
                            feasilist.Add(currentobj);
                        }
                        else
                        {
                            optlist.Add(currentobj);
                            Nooptsols++;
                        }
                    }
                }

                Program.g_comparativeIndicators.WriteLine($"{Nofeasols}, {Nooptsols}, {optobj}, {bestbound}, {relativeGap}, {totalCPUtime}");
                Program.g_comparativeIndicators.WriteLine();

                Program.g_comparativeIndicators.Write("feasiblesol" + ",");

                for (int i = 0; i < Nofeasols; i++)
                {
                    Program.g_comparativeIndicators.Write($"{feasilist[i]}" + ",");
                }
                Program.g_comparativeIndicators.WriteLine();

                Program.g_comparativeIndicators.Write("Optimalsol" + ",");

                for (int i = 0; i < Nooptsols; i++)
                {
                    Program.g_comparativeIndicators.Write($"{optlist[i]}" + ",");
                }

                Program.g_comparativeIndicators.WriteLine();

                Program.g_comparativeIndicators.WriteLine("out of memory error occurs:,{0}", e.Message);

                Program.SolutionIteration.WriteLine();

                int NumofNodesExplored = bdmaster.model.Nnodes;

                Program.g_comparativeIndicators.WriteLine();
                Program.g_comparativeIndicators.WriteLine("TreeDepth, # NodesExplored, #Iter, Obj,bestbound, Gap, CPU time");
                Program.g_comparativeIndicators.WriteLine($"{bdmaster.TreeDepth},{NumofNodesExplored},{bdmaster.NumofIter}, {optobj},{bestbound}, {relativeGap}, {totalCPUtime}");

                Cplex.Status solvests = bdmaster.model.GetStatus();

                Program.SolutionIteration.WriteLine($"Total B&C time,{totalCPUtime}");
                Program.SolutionIteration.WriteLine($"Optimal obj,{optobj}");
                Program.SolutionIteration.WriteLine($"Best Linear bound,{bestbound}");

                Console.WriteLine($"Optimal master objective cost: {bdmaster.model.GetObjValue()}");
                Console.WriteLine($"Best lower bound: {bdmaster.model.GetBestObjValue()}");

                bdmaster.model.End();

            }
            solution.output_solution(y_solution, omegaval); solution.outputSubproSolution(y_solution);
            return bdmaster;

            #endregion
        }         
    }
    class CandCG
    {
        Data data = new Data();
        Solution solution = new Solution();

        CCGSub subproblem;
        CCGMaster masterproblem = new CCGMaster();
        DualCCGSub dualsubproblem = new DualCCGSub();
        int iter_l = 0;

        public void CCGwBAC()
        {
            #region solution process
            
            Program.g_CCGiteration.WriteLine("IterNumber, UB, LB, etaval, worstscenariocost, suproblemtime, mastertime");
            Program.g_CPLEXResults.WriteLine("# feasol, # optsol, best integer, best bound, Gap, CPU time");

            List<Path> feedbacktomaster = new List<Path>();
            double eta_value = 0, mastercost = 0, UB = float.MaxValue, LB = float.MinValue;

            Stopwatch subpro = new Stopwatch();
            Stopwatch masterpro = new Stopwatch();
            Stopwatch totalProcedure = new Stopwatch();
            double mastertime = 0;
            double subtime = 0;

            double mastertime_iter = 0;
            double subtime_iter = 0;

            totalProcedure.Start();
            double availabletime = data.TL;
            int cutting_flag = 0;
            int[] y_solution = new int[data.DCSize];
            subproblem = new CCGSub();
            subproblem.GenCCGSubproblem(y_solution, new int[data.DCSize]);
            
            subproblem.model.SetParam(Cplex.Param.Threads, 1);

            masterproblem.InitializeCCGMaster();
            masterproblem.model.SetParam(Cplex.Param.Threads, 1);
            masterproblem.model.SetParam(Cplex.Param.TimeLimit, availabletime);

            CACGBD colBDGen = new CACGBD();
            
            int number_of_opens = 0;
            
            bool feasible_flag = false;
            double worst_scenario_cost = float.MinValue;

            subpro = new Stopwatch();
            subpro.Start();
            double tickcountstart = Environment.TickCount;
            double totalrunningtime = 0;
            int[] worst_u_j = new int[data.DCSize];

            subpro.Stop();

            subtime_iter = subpro.ElapsedMilliseconds / 1000;
            subtime += subtime_iter;

            Program.g_CCGiteration.WriteLine($"{iter_l},{UB}, {LB},{eta_value}, {eta_value}, {worst_scenario_cost}, {subtime_iter}, {mastertime_iter}");
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

            List<int[]> cuttingflaglist = new List<int[]>();
            List<int[]> subproblemscenarios = new List<int[]>();
            cuttingflaglist.Add(worst_u_j);

            double eps = 0.00001;
            iter_l = -1;

            while (UB - LB > eps * UB && totalrunningtime <= data.TL)
            {
                iter_l++;
                masterpro = new Stopwatch();
                masterpro.Start();

                if (cuttingflaglist.Count != 0)
                {
                    cutting_flag = 1;
                    for (int s = 0; s < cuttingflaglist.Count; s++)
                    {
                        int[] candidateuj = cuttingflaglist[s];

                        masterproblem.u_l_j.Add(candidateuj);
                        masterproblem.GenCCGMasterproblem(cutting_flag);
                        masterproblem.model.SetParam(Cplex.Param.TimeLimit, availabletime);
                    }
                    
                }
                else
                {
                    if(iter_l != 0 && !feasible_flag)
                    {
                        cutting_flag = 0;

                        masterproblem.u_l_j.Add(worst_u_j);
                        masterproblem.GenCCGMasterproblem(cutting_flag);
                        masterproblem.model.SetParam(Cplex.Param.TimeLimit, availabletime);
                    }                    
                }

                masterproblem.model.Solve();

                LB = masterproblem.model.ObjValue;
                mastercost = masterproblem.model.ObjValue;
                eta_value = masterproblem.model.GetValue(masterproblem.omega);
                
                
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
                
                masterpro.Stop();

                mastertime_iter = masterpro.ElapsedMilliseconds/1000;
                mastertime += mastertime_iter;

                feasible_flag = false;
                cuttingflaglist = new List<int[]>();

                subpro = new Stopwatch();
                subpro.Start();
                worst_u_j = new int[data.DCSize];

                
                SCNRMaster submaster = colBDGen.BranchandCutForWorstCsenario(y_solution);
                Cplex.Status cur_status;
                if (!submaster.feasiblestatus)
                {
                    worst_u_j = submaster.bestFeasibleSolution.ToArray();

                    CCGSub ccgsub = new CCGSub();

                    ccgsub.ArcFLowModel(y_solution, worst_u_j);
                    ccgsub.model.Solve();
                    double objcost = ccgsub.model.GetObjValue();
                }
                else
                {
                    cur_status = submaster.model.GetStatus();
                    
                    if (cur_status == Cplex.Status.Infeasible)
                    {
                        Console.WriteLine("Please check!!! No feasible solutions available!!!!");
                        totalrunningtime = data.TL + 1; ;
                    }
                    
                    for (int l = 0; l < data.DCSize; l++)
                    {
                        var u_val = submaster.model.GetValue(submaster.u_j[l]);

                        if (u_val > 0.5)
                        {
                            worst_u_j[l] = 1;
                        }
                    }
                    
                    subproblem.ResetOBJaRHS_XI(y_solution, worst_u_j);
                    
                    
                    subproblem.model.Solve();
                    Cplex.Status sub_status = subproblem.model.GetStatus();

                    if (sub_status == Cplex.Status.Optimal || sub_status == Cplex.Status.Feasible)
                    {
                        feasible_flag = true;
                        worst_scenario_cost = subproblem.model.GetObjValue();

                        
                        UB = Math.Min(UB, worst_scenario_cost + mastercost - eta_value);
                        
                        cuttingflaglist.Add(worst_u_j);
                                               

                        if (data.multiCutStrategy)
                        {
                            for (int s = 0; s < colBDGen.xI_sols_pool.ToList().Count; s++)
                            {
                                subproblem.ResetOBJaRHS_XI(y_solution, colBDGen.xI_sols_pool.ToList()[s]);
                                bool solved = subproblem.model.Solve();

                                sub_status = subproblem.model.GetStatus();
                                if (sub_status == Cplex.Status.Optimal || sub_status == Cplex.Status.Feasible)
                                {
                                    double ojbval = subproblem.model.GetObjValue();
                                    
                                    cuttingflaglist.Add(colBDGen.xI_sols_pool.ToList()[s]);
                                }
                                
                            }
                            
                        }
                    }
                    
                }
                
                subpro.Stop();

                subtime_iter = subpro.ElapsedMilliseconds/1000;
                subtime += subtime_iter;
                Program.g_CCGiteration.WriteLine($"{iter_l},{UB}, {LB},{eta_value}, {worst_scenario_cost}, {subtime_iter}, {mastertime_iter}");

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

                totalrunningtime = (Environment.TickCount - tickcountstart) / 1000;
                availabletime = data.TL - totalrunningtime;
            }

            totalProcedure.Stop();

            double optobj = UB;
            double bestbound = LB;
            double totalCPUtime = totalProcedure.ElapsedMilliseconds / 1000;
            double relativeGap = (UB - LB) / UB * 100;

            int totalSolutionCount = masterproblem.model.GetSolnPoolNsolns();
            int Nofeasols = 0;
            int Nooptsols = 0;

            List<double> optlist = new List<double>();
            List<double> feasilist = new List<double>();

            if (totalSolutionCount >= 1)
            {
                for (int n = 0; n < totalSolutionCount; n++)
                {
                    int[] y_val = new int[data.DCSize];
                    double[] z_val = new double[data.DCSize];
                    double eta_val = 0;

                    for (int j = 0; j < data.DCSize; j++)
                    {
                        if (masterproblem.model.GetValue(masterproblem.y_j[j], n) > 0.5)
                        {
                            y_val[j] = 1;
                        }
                        eta_val = masterproblem.model.GetValue(masterproblem.omega, n);
                    }
                    double currentobj = 0;
                    for (int j = 0; j < data.DCSize; j++)
                    {
                        currentobj += data.f_j[j] * y_val[j];
                        currentobj += data.c_j[j] * z_val[j];
                    }
                    currentobj += eta_val;

                    if (currentobj > optobj)
                    {
                        Nofeasols++;
                        feasilist.Add(currentobj);
                    }
                    else
                    {
                        optlist.Add(currentobj);
                        Nooptsols++;
                    }

                }
            }

            Program.g_CPLEXResults.WriteLine($"{Nofeasols}, {Nooptsols}, {optobj}, {bestbound}, {relativeGap}, {totalCPUtime}");
            Program.g_CPLEXResults.WriteLine();

            Program.g_CPLEXResults.Write("feasiblesol" + ",");

            for (int i = 0; i < Nofeasols; i++)
            {
                Program.g_CPLEXResults.Write($"{feasilist[i]}" + ",");
            }
            Program.g_CPLEXResults.WriteLine();

            Program.g_CPLEXResults.Write("Optimalsol" + ",");

            for (int i = 0; i < Nooptsols; i++)
            {
                Program.g_CPLEXResults.Write($"{optlist[i]}" + ",");
            }

            Program.g_CPLEXResults.WriteLine();
            Program.g_CPLEXResults.WriteLine(" #Iter, mastertime, subtime,  Obj, bestbound, Gap, CPU time");

            if (UB - LB > eps * UB)
            {
                relativeGap = (UB - LB) / UB;
                Program.g_CPLEXResults.WriteLine($"{iter_l + 1},{mastertime},{subtime}, {UB},{LB}, {relativeGap}, {totalCPUtime}");
            }
            else
            {
                Program.g_CPLEXResults.WriteLine($"{iter_l + 1},{mastertime},{subtime}, {optobj}, {bestbound}, {relativeGap}, {totalCPUtime}");
            }

            Program.g_CCGiteration.WriteLine();
            Program.g_CCGiteration.WriteLine($"Totalmastertime,{mastertime}");
            Program.g_CCGiteration.WriteLine($"Totalsubproblemtime,{subtime}");

            solution.write_solution(y_solution, eta_value);
            solution.write_subproblemsolution(subproblem);
            Console.WriteLine("UB:{0}, LB:{1}", UB, LB);
            #endregion
        }                
    }
    class Simulation
    {
        Data data = new Data();
        Solution solution = new Solution();
        
        public void TwoStgEva()
        {
            #region
            

            int[] y_sol = new int[data.DCSize];
            double[] z_sol = new double[data.DCSize];            
            double[] x_sol = new double[data.pathSize];
            
            CACGBD colBDGen = new CACGBD();

            BendersMaster Bendersmaster = colBDGen.BranchandCut();
            double omegaval = Bendersmaster.model.GetValue(Bendersmaster.omega);
            for (int j = 0; j < data.DCSize; j++)
            {
                var val = Bendersmaster.model.GetValue(Bendersmaster.y_j[j]);
                if (val > 0.5)
                {
                    y_sol[j] = 1;
                }
            }

            solution.output_solution(y_sol, omegaval); solution.outputSubproSolution(y_sol);
            Console.WriteLine(Bendersmaster.model.GetObjValue()); Console.WriteLine(Bendersmaster.model.GetBestObjValue());
            CCGSub subpro = new CCGSub();
            
            subpro.GenCCGSubproblem(y_sol, new int[data.DCSize]);
            int realcountperdistr = data.realizationlist.Count/ data.distributionSize;
            double[,] probrealizations = new double[data.distributionSize, realcountperdistr];
            int[] infeasiblenumarr = new int[data.distributionSize];
            double[] totaldemand = new double[data.distributionSize];
            double[] Underfulfilldemand = new double[data.distributionSize];
            List<double>[] FirststgAvgCost = new List<double>[data.distributionSize];

            for (int a = 0; a < data.distributionSize; a++)
            {
                FirststgAvgCost[a] = new List<double>();
            }
            int ind = 0;
            for (int d = 0; d < data.realizationlist.Count; d++)
            {
                Realization real = data.realizationlist[d];

                double[] demand = new double[data.DCSize];

                int realInd = real.sceseqno - real.distributionSeqno * realcountperdistr;                
                demand = data.h_i;
                totaldemand[real.distributionSeqno] += demand.Sum();

                z_sol = new double[data.DCSize];
                x_sol = new double[data.pathSize];
                
                subpro.ResetOBJaRHS_XI(y_sol, real.scenario);

                if (subpro.model.Solve())
                {
                    Cplex.Status solvestatus = subpro.model.GetStatus();
                    for (int j = 0; j < data.DCSize; j++)
                    {
                        z_sol[j] = subpro.model.GetValue(subpro.z_j[j]);

                        for (int r = 0; r < data.pathSize; r++)
                        {
                            x_sol[r] = subpro.model.GetValue(subpro.x_r[r]);
                        }
                    }

                    double evacost = Evaluation(y_sol, z_sol, x_sol, demand);

                    probrealizations[real.distributionSeqno, realInd] = evacost;
                }
                else
                {
                    double evacost = -10000;
                    probrealizations[real.distributionSeqno, realInd] = evacost;
                    
                    Underfulfilldemand[real.distributionSeqno] += demand.Sum();
                    infeasiblenumarr[real.distributionSeqno]++;
                    double totalcost = 0;
                    for (int j = 0; j < data.DCSize; j++)
                    {
                        totalcost += data.f_j[j] * y_sol[j];                        
                    }

                    FirststgAvgCost[real.distributionSeqno].Add(totalcost);
                }
            }

            Program.TwoStgSimuResults.WriteLine("0.1-simu, 0.2-simu, 0.3-simu");

            List<double>[] probbasedlist = new List<double>[data.distributionSize];

            for (int i = 0; i < data.distributionSize; i++)
            {
                probbasedlist[i] = new List<double>();
            }


            for (int r = 0; r < probrealizations.GetLength(1); r++)
            {
                for (int gamma = 0; gamma < data.distributionSize; gamma++)
                {
                    if (probrealizations[gamma, r] != -10000)
                    {
                        probbasedlist[gamma].Add(probrealizations[gamma, r]);
                    }
                    Program.TwoStgSimuResults.Write($"{probrealizations[gamma, r]},");
                }
                Program.TwoStgSimuResults.WriteLine();
            }

            Program.TwoStgSimuResults.WriteLine();
            Program.TwoStgSimuResults.WriteLine();

            int[] indexarr = new int[data.distributionSize];
            for (int i = 0; i < probbasedlist.Length; i++)
            {
                probbasedlist[i].Sort();
                indexarr[i] = probbasedlist[i].Count - (int)Math.Ceiling(probbasedlist[i].Count * 0.05);
            }
            
            for (int gamma = 0; gamma < data.distributionSize; gamma++)
            {
                Program.TwoStgSimuResults.Write($"{infeasiblenumarr[gamma]},");
            }
            Program.TwoStgSimuResults.WriteLine("NumInfeasi");
            
            for (int gamma = 0; gamma < data.distributionSize; gamma++)
            {
                Program.TwoStgSimuResults.Write($"{Underfulfilldemand[gamma] / totaldemand[gamma]},");
            }
            Program.TwoStgSimuResults.WriteLine("Underfulfill");

            
            for (int gamma = 0; gamma < data.distributionSize; gamma++)
            {
                if (FirststgAvgCost[gamma].Count == 0)
                {
                    Program.TwoStgSimuResults.Write("0,");
                }
                else
                {
                    Program.TwoStgSimuResults.Write($"{FirststgAvgCost[gamma].Average()},");
                }
            }
            Program.TwoStgSimuResults.WriteLine("FirstStgCost");

            
            for (int gamma = 0; gamma < data.distributionSize; gamma++)
            {
                Program.TwoStgSimuResults.Write($"{probbasedlist[gamma].Average()},");
            }

            Program.TwoStgSimuResults.WriteLine("Avg");

            
            for (int gamma = 0; gamma < data.distributionSize; gamma++)
            {
                Program.TwoStgSimuResults.Write($"{Math.Sqrt(probbasedlist[gamma].Select(val => Math.Pow(val - probbasedlist[gamma].Average(), 2)).Sum() / (probbasedlist[gamma].Count - 1))},");
            }

            Program.TwoStgSimuResults.WriteLine("Std");
            
            for (int gamma = 0; gamma < data.distributionSize; gamma++)
            {
                Program.TwoStgSimuResults.Write($"{probbasedlist[gamma][indexarr[gamma]]},");
            }
            Program.TwoStgSimuResults.WriteLine("VaR95");

            
            for (int i = 0; i < probbasedlist.Length; i++)
            {
                double totalcost = 0;
                double counter = 0;
                for (int j = indexarr[i]; j < probbasedlist[i].Count; j++)
                {
                    totalcost += probbasedlist[i][j];
                    counter++;
                }

                Program.TwoStgSimuResults.Write($"{totalcost / counter},");
            }
            Program.TwoStgSimuResults.Write($"CVaR95");
            #endregion
        }
        public double Evaluation(int[] y_sol, double[] z_sol, double[] x_sol, double[] demand)
        {
            #region
            double totalcost = 0;

            
            for (int j = 0; j < data.DCSize; j++)
            {
                totalcost += data.f_j[j] * y_sol[j];
                
            }            
            double[] u_sol = new double[data.DCSize];

            
            for (int j = 0; j < data.DCSize; j++)
            {
                totalcost += data.c_j[j] * z_sol[j];
                
            }
            for (int r = 0; r < data.pathSize; r++)
            {
                totalcost += data.d_r[r] * demand[data.e_r[r]] * x_sol[r];
            }            

            return totalcost;

            #endregion
        }
        public double EvaluationNFM(int[] y_sol, double[] z_sol, double[,] f_sol, double[] demand)
        {
            #region
            double totalcost = 0;

            
            for (int j = 0; j < data.DCSize; j++)
            {
                totalcost += data.f_j[j] * y_sol[j];

            }
            double[] u_sol = new double[data.DCSize];

            
            for (int j = 0; j < data.DCSize; j++)
            {
                totalcost += data.c_j[j] * z_sol[j];

            }
            for (int j = 0; j < data.DCSize; j++)
            {
                for (int r = 0; r < data.linkSize; r++)
                {
                    totalcost += data.linklist[r].LinkCost * f_sol[j,r];
                }
            }
            

            return totalcost;

            #endregion
        }
        public void NetworkFlow()
        {            
            #region
            

            int[] y_sol = new int[data.DCSize];
            double[] z_sol = new double[data.DCSize];
            double[,] f_sol = new double[data.DCSize,data.linkSize];

            CACGBD colBDGen = new CACGBD();

            BendersMaster Bendersmaster = colBDGen.BranchandCut();
            double omegaval = Bendersmaster.model.GetValue(Bendersmaster.omega);
            for (int j = 0; j < data.DCSize; j++)
            {
                var val = Bendersmaster.model.GetValue(Bendersmaster.y_j[j]);
                if (val > 0.5)
                {
                    y_sol[j] = 1;
                }
            }

            solution.output_solution(y_sol, omegaval); solution.outputSubproSolution(y_sol);
            Console.WriteLine(Bendersmaster.model.GetObjValue()); Console.WriteLine(Bendersmaster.model.GetBestObjValue());
            
            CCGSub subpro = new CCGSub();

            subpro.ArcFLowModel(y_sol, new int[data.DCSize]);
            int realcountperdistr = data.realizationlist.Count / data.distributionSize;
            double[,] probrealizations = new double[data.distributionSize, realcountperdistr];
            int[] infeasiblenumarr = new int[data.distributionSize];
            double[] totaldemand = new double[data.distributionSize];
            double[] Underfulfilldemand = new double[data.distributionSize];
            List<double>[] FirststgAvgCost = new List<double>[data.distributionSize];

            for (int a = 0; a < data.distributionSize; a++)
            {
                FirststgAvgCost[a] = new List<double>();
            }
            int ind = 0;

            for (int d = 0; d < data.realizationlist.Count; d++)
            {
                Realization real = data.realizationlist[d];

                double[] demand = new double[data.DCSize];

                int realInd = real.sceseqno - real.distributionSeqno * realcountperdistr;
                
                demand = data.h_i;
                totaldemand[real.distributionSeqno] += demand.Sum();

                z_sol = new double[data.DCSize];
                f_sol = new double[data.DCSize, data.linkSize];

                subpro = new CCGSub();
                subpro.ArcFLowModel(y_sol, real.scenario);

                CCGSub pathsubpro = new CCGSub();
                pathsubpro.GenCCGSubproblem(y_sol, real.scenario);

                if (subpro.model.Solve())
                {
                    Cplex.Status solvestatus = subpro.model.GetStatus();
                    for (int j = 0; j < data.DCSize; j++)
                    {
                        z_sol[j] = subpro.model.GetValue(subpro.z_j[j]);

                        for (int r = 0; r < data.linkSize; r++)
                        {
                            f_sol[j,r] = subpro.model.GetValue(subpro.f[j,r]);
                        }
                    }
                    
                    
                    

                    
                    

                    double evacost = EvaluationNFM(y_sol, z_sol, f_sol, demand);
                    
                    probrealizations[real.distributionSeqno, realInd] = evacost;
                }
                else
                {
                    double evacost = -10000;
                    probrealizations[real.distributionSeqno, realInd] = evacost;
                    
                    Underfulfilldemand[real.distributionSeqno] += demand.Sum();
                    infeasiblenumarr[real.distributionSeqno]++;
                    double totalcost = 0;
                    for (int j = 0; j < data.DCSize; j++)
                    {
                        totalcost += data.f_j[j] * y_sol[j];
                    }

                    FirststgAvgCost[real.distributionSeqno].Add(totalcost);
                }
            }

            Program.TwoStgSimuResults.WriteLine("0.1-simu, 0.2-simu, 0.3-simu");

            List<double>[] probbasedlist = new List<double>[data.distributionSize];

            for (int i = 0; i < data.distributionSize; i++)
            {
                probbasedlist[i] = new List<double>();
            }


            for (int r = 0; r < probrealizations.GetLength(1); r++)
            {
                for (int gamma = 0; gamma < data.distributionSize; gamma++)
                {
                    if (probrealizations[gamma, r] != -10000)
                    {
                        probbasedlist[gamma].Add(probrealizations[gamma, r]);
                    }
                    Program.TwoStgSimuResults.Write($"{probrealizations[gamma, r]},");
                }
                Program.TwoStgSimuResults.WriteLine();
            }

            Program.TwoStgSimuResults.WriteLine();
            Program.TwoStgSimuResults.WriteLine();

            int[] indexarr = new int[data.distributionSize];
            for (int i = 0; i < probbasedlist.Length; i++)
            {
                probbasedlist[i].Sort();
                indexarr[i] = probbasedlist[i].Count - (int)Math.Ceiling(probbasedlist[i].Count * 0.05);
            }
            
            for (int gamma = 0; gamma < data.distributionSize; gamma++)
            {
                Program.TwoStgSimuResults.Write($"{infeasiblenumarr[gamma]},");
            }
            Program.TwoStgSimuResults.WriteLine("NumInfeasi");
            
            for (int gamma = 0; gamma < data.distributionSize; gamma++)
            {
                Program.TwoStgSimuResults.Write($"{Underfulfilldemand[gamma] / totaldemand[gamma]},");
            }
            Program.TwoStgSimuResults.WriteLine("Underfulfill");

            
            for (int gamma = 0; gamma < data.distributionSize; gamma++)
            {
                if (FirststgAvgCost[gamma].Count == 0)
                {
                    Program.TwoStgSimuResults.Write("0,");
                }
                else
                {
                    Program.TwoStgSimuResults.Write($"{FirststgAvgCost[gamma].Average()},");
                }
            }
            Program.TwoStgSimuResults.WriteLine("FirstStgCost");

            
            for (int gamma = 0; gamma < data.distributionSize; gamma++)
            {
                Program.TwoStgSimuResults.Write($"{probbasedlist[gamma].Average()},");
            }

            Program.TwoStgSimuResults.WriteLine("Avg");

            
            for (int gamma = 0; gamma < data.distributionSize; gamma++)
            {
                Program.TwoStgSimuResults.Write($"{Math.Sqrt(probbasedlist[gamma].Select(val => Math.Pow(val - probbasedlist[gamma].Average(), 2)).Sum() / (probbasedlist[gamma].Count - 1))},");
            }

            Program.TwoStgSimuResults.WriteLine("Std");
            
            for (int gamma = 0; gamma < data.distributionSize; gamma++)
            {
                Program.TwoStgSimuResults.Write($"{probbasedlist[gamma][indexarr[gamma]]},");
            }
            Program.TwoStgSimuResults.WriteLine("VaR95");

            
            for (int i = 0; i < probbasedlist.Length; i++)
            {
                double totalcost = 0;
                double counter = 0;
                for (int j = indexarr[i]; j < probbasedlist[i].Count; j++)
                {
                    totalcost += probbasedlist[i][j];
                    counter++;
                }

                Program.TwoStgSimuResults.Write($"{totalcost / counter},");
            }
            Program.TwoStgSimuResults.Write($"CVaR95");
            #endregion
        }
    }
}
