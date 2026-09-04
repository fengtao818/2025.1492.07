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
        public double righthandsideU;

        public double[] x_r_val;        
        public double[] v_i_val;
        public double[] delta_r_val;
        public double[] w_j_val;

        public double[,] MXi_val;
        public DataStructure(Data data)
        {
            x_r_val = new double[data.pathSize];
            v_i_val = new double[data.nodeSize];
            MXi_val = new double[data.nodeSize, data.DCSize];
            delta_r_val = new double[data.pathSize];
            w_j_val = new double[data.DCSize];
        }
        #endregion
    }
    internal class BendersLazyConsCallback : Cplex.LazyConstraintCallback
    {
        internal BendersMaster Bendmaster;
        internal Data data;
        internal int count;
        //internal function
        internal CCGSub ccgsub;
        internal DualCCGSub dccgsub;
        internal BendersDualSub bddsub;
        bool rootnodefathom;
        //internal BendersDualSub bddsub;
        /// <summary>
        /// initialize BendersLazyConsCallback
        /// </summary>
        /// <param name="BDmaster"></param master problem with y and z>
        /// <param name="dt"></param the set of initial benders cuts>
        public BendersLazyConsCallback(BendersMaster BDmaster, Data dt)
        {
            Bendmaster = BDmaster;
            data = dt;
            count = 0;
            //generate a new CCG subproblem
            ccgsub = new CCGSub();
            ccgsub.GenCCGSubproblem(new int[data.DCSize], new int[data.DCSize]);

            bddsub = new BendersDualSub();
            bddsub.GendeltaUBsubDual(new int[data.DCSize], new int[data.DCSize]);
            bddsub.model.SetParam(Cplex.Param.Preprocessing.Presolve, false);

            dccgsub = new DualCCGSub();
            dccgsub.GenDualofCCGsub(new int[data.DCSize], new double[data.pathSize]);
            dccgsub.model.SetParam(Cplex.Param.MIP.Pool.Capacity, data.maxMultiCuts); // 解池容量
            dccgsub.model.SetParam(Cplex.Param.MIP.Pool.Replace, 2);
            dccgsub.model.SetParam(Cplex.Param.MIP.Pool.Intensity, 0);

            rootnodefathom = true;
        }
        public override void Main()
        {
            #region BranchandCheck

            double[] y_sol = new double[data.DCSize];
            double[] z_solution = new double[data.DCSize];
            int[] worst_u_j = new int[data.DCSize];
            List<int[]> multicutslist = new List<int[]>();

            Stopwatch subproblemtime = new Stopwatch();
            Stopwatch multicutstime = new Stopwatch();
            Stopwatch paretotime = new Stopwatch();
            Stopwatch findw0rstscenario = new Stopwatch();

            double subproblemtimerecord = 0;
            double multicutstimerecord = 0;
            double paretotimerecord = 0;
            double findw0rstscenariorecord = 0;

            int numberofopens = 0;
            for (int j = 0; j < data.DCSize; j++)
            {
                //worst_xi_j[j] = 0;
                y_sol[j] = GetValue(Bendmaster.y_j[j]);
                if (y_sol[j] > 0.5)
                {
                    numberofopens++;
                }                
            }
            count++;
            int[] y_solution = new int[data.DCSize];
            for (int j = 0; j < data.DCSize; j++)
            {
                if (y_sol[j] > 0.5)
                {
                    y_solution[j] = 1;
                }
            }
            if (numberofopens == 0)
            {
                for (int j = 0; j < data.DCSize; j++)
                {
                    y_solution[j] = 1 - worst_u_j[j];
                }
            }

            CCGSub ccgsubproblem = new CCGSub();
            List<Cplex> subproblemlist = new List<Cplex>();
            List<int[]> scenarios = new List<int[]>();
            List<double> costlist = new List<double>();

            object lockObject = new object();//protect the file
            double worst_scenario_cost = float.MinValue;
            
            //stg 1: getting worst scenario by enumeratedly solving subproblem 
            if (!data.solvingPrimalSub)
            {
                findw0rstscenario.Start();
                //this strategy needs to enumerate all possible scenarios which is possibly time-consuming
                List<int> open_DCent_No_list = new List<int>();
                for (int j = 0; j < data.DCSize; j++)
                {
                    if (y_solution[j] == 1)
                    {
                        open_DCent_No_list.Add(j);
                    }
                }

                // All scenarios
                List<List<int>> all_scenarios = new List<List<int>>();
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
                double[] subcostarr = new double[all_scenarios.Count];
                for (int l = 0; l < all_scenarios.Count; l++)
                {
                    List<int> temp_solution = all_scenarios[l].ToList();

                    int[] Dcent_state = new int[data.DCSize];

                    for (int j = 0; j < data.DCSize; j++)
                    {
                        Dcent_state[j] = 0;//initialization;                        
                    }

                    for (int j = 0; j < temp_solution.Count; j++)
                    {
                        int dcentno = temp_solution[j];
                        Dcent_state[dcentno] = y_solution[dcentno];//diruption                                                                   
                    }
                    ccgsubproblem = new CCGSub();
                    ccgsubproblem.GenCCGSubproblem(y_solution, Dcent_state);
                    ccgsubproblem.model.ExportModel("ccgsubproblem.lp");
                    //collect subproblem                  
                    subproblemlist.Add(ccgsubproblem.model);
                    //collect scenarios
                    scenarios.Add(Dcent_state);
                }

                bool infeasibleFound = false;  // 共享标志
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
                        lock (lockObject)
                        {
                            infeasibleFound = true;
                        }
                        state.Stop();
                    }

                });
                if (!infeasibleFound)
                {
                    for (int s = 0; s < subproblemlist.Count; s++)
                    {
                        double objcost = subproblemlist[s].GetObjValue();
                        costlist.Add(objcost);
                        subcostarr[s] = objcost;
                        if (objcost > worst_scenario_cost)
                        {
                            worst_scenario_cost = objcost;
                            worst_u_j = scenarios[s].ToArray();                            
                        }
                        subproblemlist[s].End();
                    }
                }

                findw0rstscenario.Stop(); findw0rstscenariorecord = findw0rstscenario.ElapsedMilliseconds;
            }
            
            //stg 2: getting worst scenario by branching and cut
            if (data.BACforworstscenario)
            {
                findw0rstscenario.Start();
                Dictionary<int, double> gvval = new Dictionary<int, double>();
                CACGBD colBDGen = new CACGBD();
                SCNRMaster submaster = colBDGen.BranchandCutForWorstCsenario(y_solution);

                // Find the worst scenario by BAC in the context of the current solution
                Cplex.Status cur_status = submaster.model.GetStatus();
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

                findw0rstscenario.Stop(); findw0rstscenariorecord = findw0rstscenario.ElapsedMilliseconds;               
            }
            else
            {
                findw0rstscenario.Start();

                // Find the worst scenario by MILP in the context of the current solution
                
                //determine the upper bound of delta in case of ysolution
                bddsub.ResetBDDSubObj(y_solution, y_solution);
                bddsub.model.Solve();

                double[] deltaval = new double[data.pathSize];
                deltaval = bddsub.model.GetValues(bddsub.delta_r);

                //update latest y and delta
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
            //Step: Solving ccg subproblem for dual
            subproblemtime.Start();

            //Reset the objective and right hand side of constraints related to y z xi solutions
            ccgsub.ResetOBJaRHS_XI(y_solution, worst_u_j);
            ccgsub.model.Solve();
            worst_scenario_cost = ccgsub.model.GetObjValue();

            //ccgsub.model.ExportModel("ccgsubtestmodel.lp");
            //generating optimality cut
            IRange bendersoptimalitycut = ccgsub.GenBendersCut(Bendmaster, worst_u_j);
            if (bendersoptimalitycut != null)
            {
                Add(bendersoptimalitycut);
                Bendmaster.cutsStore.Add(bendersoptimalitycut);
            }

            subproblemtime.Stop();
            subproblemtimerecord += subproblemtime.ElapsedMilliseconds;

            multicutstime.Start();

            // add multi cuts generated by local branching
            if (data.multiCutStrategy)
            {
                while (multicutslist.Count != 0)
                {
                    int[] scenario = multicutslist[0];
                    multicutslist.RemoveAt(0);

                    for (int r = 0; r < data.pathSize; r++)
                    {
                        ccgsub.relatedtodual_Delta_value[r].LB = y_solution[data.s_r[r]] * (scenario[data.s_r[r]] - 1);
                    }

                    ccgsub.model.Solve();
                    IRange optimalitycut = ccgsub.GenBendersCut(Bendmaster, scenario);

                    if (optimalitycut != null)
                    {
                        Add(optimalitycut);
                        Bendmaster.cutsStore.Add(optimalitycut);
                    }
                }
            }

            multicutstime.Stop();
            multicutstimerecord += multicutstime.ElapsedMilliseconds;

            paretotime.Start();

            //add pareto optimal subproblem cut
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
            
            Program.SolutionIteration.WriteLine($"{findw0rstscenariorecord / 1000}, {subproblemtimerecord / 1000.0}," +
                $"{multicutstimerecord / 1000.0},{paretotimerecord / 1000.0},{worst_scenario_cost},{GetBestObjValue()},{worst_u_j.ToString()}");

            #endregion
        }        
    } // END BendersLazyConsCallback
    internal class BLCCforCCGMaster : Cplex.LazyConstraintCallback
    {
        internal CCGMaster Bendmaster;
        internal Data data;
        internal int count;
        //internal function
        internal CCGSub ccgsub;
        
        bool rootnodefathom;
        List<int[]> multicutslist;
        //internal BendersDualSub bddsub;
        /// <summary>
        /// initialize BendersLazyConsCallback
        /// </summary>
        /// <param name="BDmaster"></param master problem with y and z>
        /// <param name="dt"></param the set of initial benders cuts>
        public BLCCforCCGMaster(CCGMaster BDmaster, Data dt, List<int[]> multiscenario)
        {
            Bendmaster = BDmaster;
            data = dt;
            multicutslist = multiscenario;
            count = 0;
            //generate a new CCG subproblem
            ccgsub = new CCGSub();
            ccgsub.GenCCGSubproblem(new int[data.DCSize], new int[data.DCSize]);

            rootnodefathom = true;
        }
        public override void Main()
        {
            #region BranchandCheck

            double[] y_sol = new double[data.DCSize];
            double[] z_solution = new double[data.DCSize];
            int[] worst_u_j = new int[data.DCSize];
            
            Stopwatch subproblemtime = new Stopwatch();
            Stopwatch multicutstime = new Stopwatch();
            Stopwatch paretotime = new Stopwatch();
            Stopwatch findw0rstscenario = new Stopwatch();

            double subproblemtimerecord = 0;
            double multicutstimerecord = 0;
            double paretotimerecord = 0;
            double findw0rstscenariorecord = 0;

            int numberofopens = 0;
            for (int j = 0; j < data.DCSize; j++)
            {
                //worst_xi_j[j] = 0;
                y_sol[j] = GetValue(Bendmaster.y_j[j]);
                if (y_sol[j] > 0.5)
                {
                    numberofopens++;
                }
            }
            count++;
            int[] y_solution = new int[data.DCSize];
            for (int j = 0; j < data.DCSize; j++)
            {
                if (y_sol[j] > 0.5)
                {
                    y_solution[j] = 1;
                }
            }
            
            multicutstime.Start();
            double worst_scenario_cost = float.MinValue;
            // add multi cuts generated by local branching
            if (data.multiCutStrategy)
            {
                for (int s = 0; s < multicutslist.Count; s++)
                {
                    int[] scenario = multicutslist[s];
                    //multicutslist.RemoveAt(0);

                    for (int r = 0; r < data.pathSize; r++)
                    {
                        ccgsub.relatedtodual_Delta_value[r].LB = y_solution[data.s_r[r]] * (scenario[data.s_r[r]] - 1);
                    }

                    ccgsub.model.Solve();

                    if (ccgsub.model.GetObjValue() > worst_scenario_cost)
                    {
                        worst_scenario_cost = ccgsub.model.GetObjValue();
                        worst_u_j = scenario.ToArray();
                    }

                    IRange optimalitycut = ccgsub.GenBendersCutforCCG(Bendmaster, scenario);

                    if (optimalitycut != null)
                    {
                        Add(optimalitycut);
                        Bendmaster.cutsStore.Add(optimalitycut);
                    }
                }                
            }

            multicutstime.Stop();
            multicutstimerecord += multicutstime.ElapsedMilliseconds;

            paretotime.Start();

            //add pareto optimal subproblem cut
            if (data.paretoCutStrategy)
            {
                BendersDualSub pareto = new BendersDualSub();
                pareto.GenBendersubDual(y_solution, worst_u_j);
                IRange paretoCut = pareto.GenParetoCutforCCGMaster(y_solution, worst_u_j, worst_scenario_cost, data.pareto_y_sol, data.pareto_u_sol, Bendmaster);

                if (paretoCut != null)
                {
                    Add(paretoCut);
                    Bendmaster.cutsStore.Add(paretoCut);
                }
                pareto.model.End(); pareto = null;
            }

            paretotime.Stop();
            paretotimerecord += paretotime.ElapsedMilliseconds;

            Program.SolutionIteration.WriteLine($"{findw0rstscenariorecord / 1000}, {subproblemtimerecord / 1000.0}," +
                $"{multicutstimerecord / 1000.0},{paretotimerecord / 1000.0},{worst_scenario_cost},{GetBestObjValue()},{worst_u_j.ToString()}");

            #endregion
        }
    } // END BendersLazyConsCallback
    /// <summary>
    /// invoked at node which finds a linear solution
    /// </summary>
    internal class BendersUserCutCallback : Cplex.UserCutCallback
    {
        internal CCGMaster Bendmaster;
        internal Data data;
        //internal function
        internal CCGSub ccgsub;
        internal int Iteration;
        public BendersUserCutCallback(CCGMaster BDmaster, Data dt)
        {
            Bendmaster = BDmaster;
            data = dt;
            Iteration = 0;
            //generate a new CCG subproblem
            ccgsub = new CCGSub();
            ccgsub.GenlinearSubproblem(new double[data.DCSize], new int[data.DCSize]);
        }
        public override void Main()
        {
            #region BranchandCheck

            if (!IsAfterCutLoop())
                return;

            double[] y_sol = new double[data.DCSize];
            double[] z_solution = new double[data.DCSize];

            for (int j = 0; j < data.DCSize; j++)
            {
                //worst_xi_j[j] = 0;
                y_sol[j] = GetValue(Bendmaster.y_j[j]);                
            }

            //Step: Solving ccg subproblem for dual

            //Reset the objective and right hand side of constraints related to y z xi solutions
            ccgsub.ResetLinearOBJaRHS_XI(y_sol, new int[data.DCSize]);
            ccgsub.model.Solve();

            IRange bendersoptimalitycut = ccgsub.GenBendersCutforCCG(Bendmaster, new int[data.DCSize]);
            if (bendersoptimalitycut != null)
            {
                AddLocal(bendersoptimalitycut);
                Bendmaster.cutsStore.Add(bendersoptimalitycut);
            }

            Iteration++;

            if (GetCurrentNodeDepth() == 0)
            {
                if (Iteration > 2000)
                {
                    Iteration = 0;
                    AbortCutLoop();
                }
            }
            else
            {
                if (Iteration > 20)
                {
                    Iteration = 0;
                    AbortCutLoop();
                }
            }
            #endregion
        }
        public bool ArraysAreEqual(int[] array1, int[] array2)
        {
            #region
            for (int j = 0; j < data.DCSize; j++)
            {
                if (array1[j] != array2[j])
                    return false;
            }
            return true;
            #endregion
        }
    } // END BendersLazyConsCallback
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
            
            //reset the objective of subproblem and the right hand side of corresponding constraints
            modifiedSub.ResetOBJaRHS_XI(y_solution, u_sol, data.big_M_r);
            bool issolved = modifiedSub.model.Solve();
            //modifiedSub.model.ExportModel("modifiedSub.lp");
            Cplex.Status solvingstatus = modifiedSub.model.GetStatus();

            if (issolved)
            {
                double[] x_r_values = new double[data.pathSize];
                bool isfeaisble = true;
                for (int r = 0; r < data.pathSize; r++)
                {
                    x_r_values[r] = modifiedSub.model.GetValue(modifiedSub.x_r[r]);
                    if (u_sol[data.s_r[r]] * x_r_values[r] >= data.epsilon)
                    {
                        Program.Solvinglog.WriteLine("Error during finding the scenario using lazyconstraint callback: " + x_r_values[r]);
                        Scenariomaster.worstScenarioSolution = u_sol.ToArray();
                        Scenariomaster.feasiblestatus = false;
                        
                        isfeaisble = false;
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
                        for (int j = 0; j < data.DCSize; j++)
                        {
                            if (Scenariomaster.bestFeasibleSolution[j] != u_sol[j])
                            {
                                sametobestsolution = false; break;
                            }
                        }

                        if (data.nogoodcuts && benders_primal_cost < Scenariomaster.bestObjVal)
                        {
                            ILinearNumExpr combinatorialCut = Scenariomaster.model.LinearNumExpr();

                            double cardi = 0;
                            for (int j = 0; j < data.DCSize; j++)
                            {
                                if (u_sol[j] == 1)
                                {
                                    combinatorialCut.AddTerm(-1, Scenariomaster.u_j[j]);
                                    cardi++;
                                }
                                else
                                {
                                    combinatorialCut.AddTerm(1, Scenariomaster.u_j[j]);
                                }
                            }
                            Add(Scenariomaster.model.Ge(combinatorialCut, 1 - cardi));

                        }
                    }

                    ILinearNumExpr extrem_point_exp = Scenariomaster.model.LinearNumExpr();

                    for (int r = 0; r < data.pathSize; r++)
                    {                        
                        extrem_point_exp.AddTerm(-data.big_M_r[r] * x_r_values[r], Scenariomaster.u_j[data.s_r[r]]);
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
                    for (int j = 0; j < data.DCSize; j++)
                    {
                        if (u_sol[j] == 1)
                        {
                            combinatorialCut.AddTerm(-1, Scenariomaster.u_j[j]);
                            cardi++;
                        }
                        else
                        {
                            combinatorialCut.AddTerm(1, Scenariomaster.u_j[j]);
                        }
                    }
                    Add(Scenariomaster.model.Ge(combinatorialCut, 1 - cardi));
                }
                if(!Scenariomaster.feasiblestatus)
                    Abort();
            }
            
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
        //used for checking the correctness of B&C
        public int[] worstScenarioSolution;
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
        /// <summary>
        /// generate MILP model for scenario master problem
        /// </summary>
        /// <param name="y_sol"></param>
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
                // set the upperbound of xi_j as the value of y_j
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
                //add terms to objective                   
                extrem_point_exp.AddTerm(1, eta);
                //add benders cut
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
        //internal INumVar[] s_i;

        public IRange[] relatedtodual_V_value;
        public IRange[] relatedtodual_W_value;
        public IRange[] relatedtodual_Delta_value;
        public IRange[] relatedtodual_Gamma_value;

        public int number_of_var;
        public int number_of_con;

        internal int[] u_sol;
        internal int[] y_sol;
        
        /// <param name="falocsolution"></param>        
        /// <param name="scenario_xi_j"></param>
        public SCNRSub(int[] y_sol, int[] u_sol, Data data)
        {
            this.y_sol = y_sol;            
            this.u_sol = u_sol;
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
            //s_i = new INumVar[data.nodeSize];
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
            //for (int i = 0; i < data.nodeSize; i++)
            //{
            //    s_i[i] = model.NumVar(0, 1, NumVarType.Float, $"s_{i}");
            //    number_of_var++;
            //}
            ILinearNumExpr singleterm = model.LinearNumExpr();
            for (int j = 0; j < data.DCSize; j++)
            {                
                singleterm.AddTerm(data.c_j[j], z_j[j]);
                
            }

            for (int r = 0; r < data.pathSize; r++)
            {                
                singleterm.AddTerm((data.d_r[r] * data.h_i[data.e_r[r]] + data.big_M_r[r] * u_sol[data.s_r[r]]), x_r[r]);
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
            if (cons_4 == 1 && data.Capacitystrategy)
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
        /// <summary>
        /// Update model at each branch or iteration without newly formulating
        /// </summary>
        /// <param name="upd_y_solution"></param-difference>        
        /// <param name="upd_u_sol"></param-difference>
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
                singleterm.AddTerm((data.d_r[r] * data.h_i[data.e_r[r]] + bigm_R[r] * upd_u_sol[data.s_r[r]]), x_r[r]);
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

        Data data = new Data();

        /// <summary>
        /// benders master problem
        /// used for the ccg master problem
        /// due to not convergent for the large instances
        /// </summary>        
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

            omega = model.NumVar(-float.MaxValue, float.MaxValue, NumVarType.Float, "omega");
            for (int j = 0; j < data.DCSize; j++)
            {
                y_j[j] = model.BoolVar($"y_{j}");
                number_of_var++; 
            }

            //only once add
            //objective function
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
            #endregion
        }
        /// <summary>
        /// add new valid cuts identified in preprocessing, i.e., solving relaxation at root node,
        /// the linear valid cuts obtained by considering stable solution
        /// </summary>
        /// <param name="coeffList"></param>
        public void AddNewValidCuts(List<DataStructure> coeffList)
        {
            #region
            int[] u_sol = new int[data.DCSize];

            for (int d = 0; d < coeffList.Count; d++)
            {
                DataStructure dts = coeffList[d];
                ILinearNumExpr extrem_point_exp = model.LinearNumExpr();

                extrem_point_exp.AddTerm(1, omega);

                //first part: y-related
                for (int r = 0; r < data.pathSize; r++)
                {
                    extrem_point_exp.AddTerm(dts.delta_r_val[r] * (1 - u_sol[data.s_r[r]]), y_j[data.s_r[r]]);
                }

                //second part: constant
                double sumofvi = 0;
                for (int i = 0; i < data.nodeSize; i++)
                {
                    sumofvi += dts.v_i_val[i];
                }
                //add Benders' cut separation
                model.AddGe(extrem_point_exp, sumofvi);
            }
            #endregion
        }        
        /// <summary>
        /// add cutoff constraints
        /// </summary>
        /// <param name="objval"></param best-known obj val>
        /// <param name="theta"></param cutoff tolerance>
        public void AddcutoffCons(double objval, double theta)
        {
            #region
            //objective function
            INumExpr numExpr = model.NumExpr();

            for (int j = 0; j < data.DCSize; j++)
            {
                numExpr = model.Sum(model.Prod(data.f_j[j], y_j[j]), numExpr);
            }
            numExpr = model.Sum(omega, numExpr);
            model.AddLe(numExpr, objval - theta, "cutoff constraint");
            #endregion
        }
        /// <summary>
        /// local branching search after rootnode
        /// </summary>
        /// <param name="incumbent"></param>
        /// <param name="radiusk"></param>
        /// <param name="branch"></param>
        public void AddlocalBranchingCons(int[] incumbent, int radiusk, int branch)
        {
            #region 

            if(branch == 0)//left
            {
                int constant = 0;
                ILinearNumExpr consexpr = model.LinearNumExpr();
                for (int j = 0; j < data.DCSize; j++)
                {                    
                    if(incumbent[j] == 1)
                    {
                        consexpr.AddTerm(-1, y_j[j]);
                        constant++;
                    }
                    else
                    {
                        consexpr.AddTerm(1, y_j[j]);
                    }
                }
                model.AddLe(consexpr, radiusk - constant, "left branch");
            }
            else//right
            {
                int constant = 0;
                ILinearNumExpr consexpr = model.LinearNumExpr();
                for (int j = 0; j < data.DCSize; j++)
                {
                    if (incumbent[j] == 1)
                    {
                        consexpr.AddTerm(-1, y_j[j]);
                        constant++;
                    }
                    else
                    {
                        consexpr.AddTerm(1, y_j[j]);
                    }
                }
                model.AddGe(consexpr, radiusk + 1 - constant, "right branch");
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
        /// <summary>
        /// ccg master problem
        /// </summary>
        /// <param name="cutting_indicator"></param>
        public void GenCCGMasterproblem(int cutting_indicator)
        {
            #region model
            //add the newly finded scenario into CCG model
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
                //add terms to objective
                model.AddGe(omega, singleterm, "optimality cut");
            }

            //constraints            
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
            if (cons_4 == 1 && data.Capacitystrategy)
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
        public void GenParetoCorePoint()
        {
            #region
            IObjective ccgmstobj = model.GetObjective();
            ccgmstobj.ClearExpr();

            ILinearNumExpr numExpr = model.LinearNumExpr();
            for (int j = 0; j < data.DCSize; j++)
            {
                numExpr.AddTerm(data.f_j[j], y_j[j]);
            }
            ccgmstobj.Expr = numExpr;


            #endregion
        }
        public void LinearCCGMaster()
        {
            #region model
            model = new Cplex();

            number_of_var = 0; number_of_con = 0;

            YCoeff_list = new List<DataStructure>();
            Temp_YCoeff_list = new List<DataStructure>();
            generateCuts = new List<IRange>();

            y_j = new INumVar[data.DCSize];
            x_l_r = new List<INumVar[]>();
            z_l_j = new List<INumVar[]>();
            u_l_j = new List<int[]>();
            
            for (int j = 0; j < data.DCSize; j++)
            {
                y_j[j] = model.NumVar(0, float.MaxValue, NumVarType.Float, $"y_{j}");
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
            int cons_1 = 1;
            if (cons_1 == 1)
            {
                ILinearNumExpr consExpr = model.LinearNumExpr();
                for (int j = 0; j < data.DCSize; j++)
                {
                    consExpr.AddTerm(1, y_j[j]);
                }
                model.AddGe(consExpr, data.max_dstroyed_DCs + 1);
            }
            int cons_2 = 1;
            if (cons_2 == 1)
            {

                for (int j = 0; j < data.DCSize; j++)
                {
                    model.AddLe(y_j[j], 1, $"upper bound_{j}");
                }

            }

            #endregion
        }
        /// <summary>
        /// benders decomposition for ccg master
        /// solving linear model (linear y) without considering \xi_j
        /// generate a set of initially valid cuts for masters
        /// </summary>
        /// <param name="xi_solution"></param>        
        public void BendersDec(double lambda, double eps_threshod, double alpha)
        {
            #region 
            double[] stab_y_sol = new double[data.DCSize];            
            //generate stable y through solving model

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

            //generate a linear model, containing linear y and other variables
            LinearCCGMaster();
            u_l_j.Add(u_sol);
            GenCCGMasterproblem(1);

            int It_counter = 0;
            //double lambda = data.Para_STB_Lambada;
            //double eps_threshod = data.Para_STB_Epsilon;
            //double alpha = data.Para_STB_Alpha;

            //The benders decomposition is based on the laster worst scenario
            bool masterfeasible = model.Solve();
            mastercost = model.GetObjValue();
            LB = mastercost;

            CCGSub ccgsub = new CCGSub();

            benders_primal_cost = -float.MaxValue;

            while (true)
            {
                double[] ast_y_sol = new double[data.DCSize];                
                double[] one_matr = new double[data.DCSize];//unit, 1

                for (int j = 0; j < data.DCSize; j++)
                {
                    //worst_xi_j[j] = 0;
                    ast_y_sol[j] = model.GetValue(y_j[j]);                    
                    one_matr[j] = 1;
                }

                // Get the current y solution
                //intermediate point

                if (It_counter >= data.Para_STB_IterLimit)
                {
                    //i.e., use primal optimal solution
                    It_counter = 0;

                    //clear
                    List<IRange> tempcutlist = new List<IRange>();
                    List<DataStructure> tempcoeff = new List<DataStructure>();
                    for (int c = 0; c < generateCuts.Count; c++)
                    {
                        double slackvalue = model.GetSlack(generateCuts[c]);
                        if (slackvalue <= 0)
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
                    //strategy 1:
                    for (int j = 0; j < data.DCSize; j++)
                    {
                        stab_y_sol[j] = alpha * stab_y_sol[j] + (1 - alpha) * ast_y_sol[j];
                        
                        y_sol[j] = lambda * ast_y_sol[j] + (1 - lambda) * stab_y_sol[j];                       
                    }
                }
                else
                {
                    //strategy 1: λxi∗ +(1− λ) ˜xi + δ(1, . . . , 1), 
                    for (int j = 0; j < data.DCSize; j++)
                    {
                        stab_y_sol[j] = (ast_y_sol[j] + stab_y_sol[j]) / 2;
                        
                        y_sol[j] = lambda * ast_y_sol[j] + (1 - lambda) * stab_y_sol[j] /*+ eps_threshod * one_matr[j]*/;                        
                    }
                }

                Cplex.Status primalstatus = null;

                ccgsub = new CCGSub();
                ccgsub.GenlinearSubproblem(y_sol, u_sol);
                ccgsub.model.ExportModel("GenDualCCGwXI.lp");
                ccgsub.model.Solve();
                benders_primal_cost = ccgsub.model.GetObjValue();
                primalstatus = ccgsub.model.GetStatus();

                if (LB > best_LP_bound)
                {
                    best_LP_bound = LB;
                    It_counter = 0;
                }
                else
                {
                    It_counter++;
                }

                if (primalstatus == Cplex.Status.Optimal)
                {
                    DataStructure dts = new DataStructure(data);

                    double[] v_i_values = new double[data.nodeSize];//>0
                    double[] delta_r_values = new double[data.pathSize];//>0
                    
                    double sumofv = 0;
                    for (int i = 0; i < data.nodeSize; i++)
                    {
                        v_i_values[i] = ccgsub.model.GetDual(ccgsub.relatedtodual_V_value[i]);
                        sumofv += v_i_values[i];
                    }

                    for (int r = 0; r < data.pathSize; r++)
                    {
                        delta_r_values[r] = Math.Abs( ccgsub.model.GetDual(ccgsub.relatedtodual_Delta_value[r]));
                    }

                    dts.delta_r_val = delta_r_values;
                    dts.v_i_val = v_i_values;
                    
                    Temp_YCoeff_list.Add(dts);

                    ILinearNumExpr extrem_point_exp = model.LinearNumExpr();

                    extrem_point_exp.AddTerm(1, omega);

                    //first part: y-related
                    for (int r = 0; r < data.pathSize; r++)
                    {
                        extrem_point_exp.AddTerm(delta_r_values[r] * (1 - u_sol[data.s_r[r]]), y_j[data.s_r[r]]);
                    }
                    generateCuts.Add(model.AddGe(extrem_point_exp, sumofv));                    
                }

                masterfeasible = model.Solve();
                mastercost = model.GetObjValue();
                LB = mastercost;

                iteration++;
            }
            #endregion
        }
        /// <summary>
        /// add new valid cuts identified in preprocessing, i.e., solving relaxation at root node,
        /// the linear valid cuts obtained by considering stable solution
        /// </summary>
        /// <param name="coeffList"></param>
        public void AddNewValidCuts(List<DataStructure> coeffList)
        {
            #region
            int[] u_sol = new int[data.DCSize];

            for (int d = 0; d < coeffList.Count; d++)
            {
                DataStructure dts = coeffList[d];
                ILinearNumExpr extrem_point_exp = model.LinearNumExpr();

                extrem_point_exp.AddTerm(1, omega);

                //first part: y-related
                for (int r = 0; r < data.pathSize; r++)
                {
                    extrem_point_exp.AddTerm(dts.delta_r_val[r] * (1 - u_sol[data.s_r[r]]), y_j[data.s_r[r]]);
                }

                //second part: constant
                double sumofvi = 0;
                for (int i = 0; i < data.nodeSize; i++)
                {
                    sumofvi += dts.v_i_val[i];
                }
                //add Benders' cut separation
                model.AddGe(extrem_point_exp, sumofvi);
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
        internal INumVar[] s_i;

        public IRange[] relatedtodual_V_value;
        public IRange[] relatedtodual_W_value;
        public IRange[] relatedtodual_Delta_value;
        public IRange[] relatedtodual_Pi_value;

        public int number_of_var;
        public int number_of_con;
        
        /// <summary>
        /// ccg subproblem 
        /// used for the valid cut for ccg master problem
        /// </summary>        
        public void GenCCGSubproblem(int[] y_sol, int[] u_sol)
        {
            #region ccg subproblem model           
            
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
            
            //add terms to objective
            model.AddMinimize(singleterm);

            //constraints         
            
            //constraints            
            int cons_1 = 1;
            if (cons_1 == 1)
            {
                for (int r = 0; r < data.pathSize; r++)
                {
                    ILinearNumExpr cons = model.LinearNumExpr();
                    cons.AddTerm(-1, x_r[r]);
                    IRange constraint = model.AddGe(cons, y_sol[data.s_r[r]] * (u_sol[data.s_r[r]] - 1), $"Transport_capacity_{r}_{data.s_r[r]}");
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
                    for (int r = 0; r < data.pathSize; r++)
                    {                        
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
            if(cons_4 == 1 && data.Capacitystrategy)
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
        /// <summary>
        /// Update model at each branch or iteration without newly formulating
        /// </summary>
        /// <param name="upd_y_solution"></param-difference>        
        /// <param name="upd_u_sol"></param-difference>
        public void ResetOBJaRHS_XI(int[] upd_y_solution, int[] upd_u_sol)
        {
            #region

            //Reset right hand side of constraints related to y z xi solutions
            for (int r = 0; r < data.pathSize; r++)
            {
                if (relatedtodual_Delta_value[r] != null)
                {
                    relatedtodual_Delta_value[r].LB = upd_y_solution[data.s_r[r]] * (upd_u_sol[data.s_r[r]] - 1);
                }
            }
            #endregion
        }
        /// <summary>
        /// generate benders cut
        /// </summary>
        /// <param name="BDM"></param master problem>
        /// <param name="xi_sol"></param worst scenario>
        /// <returns></returns>
        public IRange GenBendersCut(BendersMaster BDM, int[] u_sol)
        {
            #region
            double[] v_i_values = new double[data.nodeSize];//>0
            double[] delta_r_values = new double[data.pathSize];//>0
            double[] pi_a_values = new double[data.linkSize];

            double benders_primal_cost = model.GetObjValue();
            double sumofv = 0;
            for (int i = 0; i < data.nodeSize; i++)
            {
                v_i_values[i] = model.GetDual(relatedtodual_V_value[i]);
                sumofv += v_i_values[i];
            }

            for (int r = 0; r < data.pathSize; r++)
            {
                delta_r_values[r] = Math.Abs(model.GetDual(relatedtodual_Delta_value[r]));
            }

            if (data.Capacitystrategy)
            {
                for (int a = 0; a < data.linkSize; a++)
                {
                    pi_a_values[a] = Math.Abs(model.GetDual(relatedtodual_Pi_value[a]));

                    sumofv -= pi_a_values[a] * data.linklist[a].LinkCapacity;
                }
            }            

            ILinearNumExpr extrem_point_exp = model.LinearNumExpr();

            extrem_point_exp.AddTerm(1, BDM.omega);

            //first part: y-related
            for (int r = 0; r < data.pathSize; r++)
            {
                extrem_point_exp.AddTerm(delta_r_values[r] * (1 - u_sol[data.s_r[r]]), BDM.y_j[data.s_r[r]]);
            }
            
            IRange cut = BDM.model.Ge(extrem_point_exp, sumofv);

            return cut;
            #endregion
        }
        /// <summary>
        /// generate benders cut
        /// </summary>
        /// <param name="CCGM"></param master problem>
        /// <param name="xi_sol"></param worst scenario>
        /// <returns></returns>
        public IRange GenBendersCutforCCG(CCGMaster CCGM, int[] u_sol)
        {
            #region
            double[] v_i_values = new double[data.nodeSize];//>0
            double[] delta_r_values = new double[data.pathSize];//>0
            double[] pi_a_values = new double[data.linkSize];

            double benders_primal_cost = model.GetObjValue();
            double sumofv = 0;
            for (int i = 0; i < data.nodeSize; i++)
            {
                v_i_values[i] = model.GetDual(relatedtodual_V_value[i]);
                sumofv += v_i_values[i];
            }

            for (int r = 0; r < data.pathSize; r++)
            {
                delta_r_values[r] = Math.Abs(model.GetDual(relatedtodual_Delta_value[r]));
            }

            if (data.Capacitystrategy)
            {
                for (int a = 0; a < data.linkSize; a++)
                {
                    pi_a_values[a] = Math.Abs(model.GetDual(relatedtodual_Pi_value[a]));

                    sumofv -= pi_a_values[a] * data.linklist[a].LinkCapacity;
                }
            }

            ILinearNumExpr extrem_point_exp = model.LinearNumExpr();

            extrem_point_exp.AddTerm(1, CCGM.omega);

            //first part: y-related
            for (int r = 0; r < data.pathSize; r++)
            {
                extrem_point_exp.AddTerm(delta_r_values[r] * (1 - u_sol[data.s_r[r]]), CCGM.y_j[data.s_r[r]]);
            }

            IRange cut = CCGM.model.Ge(extrem_point_exp, sumofv);

            return cut;
            #endregion
        }
        /// <summary>
        /// linear ccg subproblem 
        /// used for the valid cut for ccg master problem
        /// </summary>
        /// <param name="y_sol"></param>        
        /// <param name="u_sol"></param>
        public void GenlinearSubproblem(double[] y_sol, int[] u_sol)
        {
            #region ccg subproblem model           

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

            //add terms to objective
            model.AddMinimize(singleterm);

            //constraints         

            //constraints            
            int cons_1 = 1;
            if (cons_1 == 1)
            {
                for (int r = 0; r < data.pathSize; r++)
                {
                    ILinearNumExpr cons = model.LinearNumExpr();
                    cons.AddTerm(-1, x_r[r]);
                    IRange constraint = model.AddGe(cons, y_sol[data.s_r[r]] * (u_sol[data.s_r[r]] - 1), $"Transport_capacity_{r}_{data.s_r[r]}");
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
                    for (int r = 0; r < data.pathSize; r++)
                    {
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
            if (cons_4 == 1 && data.Capacitystrategy)
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
        /// <summary>
        /// Update model at each branch or iteration without newly formulating
        /// </summary>
        /// <param name="upd_y_solution"></param-difference>        
        /// <param name="upd_u_sol"></param-difference>
        public void ResetLinearOBJaRHS_XI(double[] upd_y_solution, int[] upd_u_sol)
        {
            #region

            //Reset right hand side of constraints related to y z xi solutions
            for (int r = 0; r < data.pathSize; r++)
            {
                if (relatedtodual_Delta_value[r] != null)
                {
                    relatedtodual_Delta_value[r].LB = upd_y_solution[data.s_r[r]] * (upd_u_sol[data.s_r[r]] - 1);
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
        public INumVar[] u_j;
        public INumVar[] delta_r;
        public INumVar[] b_r;
        public INumVar[] pi_a;
        public IRange[] relatedtoUB_delta_value;

        public int number_of_var;
        public int number_of_con;
        
        /// <summary>
        /// dual problem of the ccg subproblem
        /// used for generating the worst scenario \xi_j
        /// note that this dual subproblem cannot provide valid cuts
        /// </summary>        
        public void GenDualofCCGsub(int[] y_sol, double[] UB_delta)
        {
            #region
            model = new Cplex();
            number_of_var = 0; number_of_con = 0;
            //model.SetParam(Cplex.Param.Preprocessing.Presolve, false);
            v_i = new INumVar[data.nodeSize];
            w_j = new INumVar[data.DCSize];
            u_j = new INumVar[data.DCSize];
            pi_a = new INumVar[data.linkSize];
            delta_r = new INumVar[data.pathSize];
            b_r = new INumVar[data.pathSize];
            relatedtoUB_delta_value = new IRange[data.pathSize];

            for (int r = 0; r < data.pathSize; r++)
            {
                delta_r[r] = model.NumVar(0, float.MaxValue, NumVarType.Float, $"delta_{r}");
                number_of_var++;
                b_r[r] = model.NumVar(0, float.MaxValue, NumVarType.Float, $"b_{r}");
                number_of_var++;
            }
            for (int i = 0; i < data.nodeSize; i++)
            {
                v_i[i] = model.NumVar(0, float.MaxValue, NumVarType.Float, $"v_{i}");
                number_of_var++;
            }
            for (int j = 0; j < data.DCSize; j++)
            {
                w_j[j] = model.NumVar(float.MinValue, 0, NumVarType.Float, $"w_{j}");
                number_of_var++;
                u_j[j] = model.NumVar(0, y_sol[j], NumVarType.Int, $"u_{j}");
                number_of_var++;
            }
            

            //first part            
            ILinearNumExpr subobj = model.LinearNumExpr();
            for (int r = 0; r < data.pathSize; r++)
            {
                subobj.AddTerm(-y_sol[data.s_r[r]], delta_r[r]);
            }
            for (int r = 0; r < data.pathSize; r++)
            {
                subobj.AddTerm(y_sol[data.s_r[r]], b_r[r]);
            }
            for (int i = 0; i < data.nodeSize; i++)
            {
                subobj.AddTerm(1, v_i[i]);
            }

            if (data.Capacitystrategy)
            {
                for (int a = 0; a < data.linkSize; a++)
                {
                    pi_a[a] = model.NumVar(0, float.MaxValue, NumVarType.Float, $"pi_{a}");
                    number_of_var++;
                }

                for (int a = 0; a < data.linkSize; a++)
                {
                    subobj.AddTerm(- data.linklist[a].LinkCapacity, pi_a[a]);
                }
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

                    consexpr.AddTerm(-1, delta_r[r]);

                    if (data.Capacitystrategy)
                    {
                        for (int l = 0; l < data.linkSize; l++)
                        {
                            if (data.pathlist[r].Contains(l))
                            {
                                consexpr.AddTerm(- data.h_i[data.e_r[r]], pi_a[l]);
                            }
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
                    IRange constraint = model.AddGe(w_j[j], -data.c_j[j], $"Z_Dual_{j}");
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

                    IRange constraint = model.AddLe(consexpr, 0, $"BigM_dual_{r}");
                    relatedtoUB_delta_value[r] = constraint;
                }                
            }
            int cons_3 = 1;
            if (cons_3 == 1)
            {
                for (int r = 0; r < data.pathSize; r++)
                {
                    model.AddLe(b_r[r], delta_r[r], $"BigM2_dual_{r}");
                    
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
        /// <summary>
        /// 
        /// </summary>
        /// <param name="upd_y_sol"></param>
        /// <param name="UB_delta"></param>
        public void ResetdualccgSubObjcons(int[] upd_y_sol, double[] upd_UB_delta)
        {
            #region
            //update the upper bound of u_j
            for (int j = 0; j < data.DCSize; j++)
            {
                u_j[j].UB = upd_y_sol[j];
            }
            //clear the objexpr in previous iteration
            IObjective originalobj = model.GetObjective();
            originalobj.ClearExpr();

            //update objective function
            ILinearNumExpr subobj = model.LinearNumExpr();
            for (int r = 0; r < data.pathSize; r++)
            {
                subobj.AddTerm(-upd_y_sol[data.s_r[r]], delta_r[r]);
            }
            for (int r = 0; r < data.pathSize; r++)
            {
                subobj.AddTerm(upd_y_sol[data.s_r[r]], b_r[r]);
            }
            for (int i = 0; i < data.nodeSize; i++)
            {
                subobj.AddTerm(1, v_i[i]);
            }
            if (data.Capacitystrategy)
            {                
                for (int a = 0; a < data.linkSize; a++)
                {
                    subobj.AddTerm(- data.linklist[a].LinkCapacity, pi_a[a]);
                }
            }
            originalobj.Expr = subobj;

            //update constraint
            for (int r = 0; r < data.pathSize; r++)
            {
                //clear the expr in previous iteration
                relatedtoUB_delta_value[r].ClearExpr();
                //write the current expr
                ILinearNumExpr consexpr = model.LinearNumExpr();
                consexpr.AddTerm(1, b_r[r]);
                consexpr.AddTerm(-upd_UB_delta[r], u_j[data.s_r[r]]);
                //update it
                relatedtoUB_delta_value[r].Expr = consexpr;
            }
            #endregion
        }
        /// <summary>
        /// local branching search after rootnode
        /// </summary>
        /// <param name="incumbent"></param>
        /// <param name="radiusk"></param should be between [data.max_dstroyed_DCs, 2 * data.max_dstroyed_DCs]>
        /// <param name="branch"></param>
        public void AddlocalBranchingCons(int[] incumbent, int radiusk, int branch)
        {
            #region 
            
            if (branch == 0)//left
            {
                int constant = 0;
                ILinearNumExpr consexpr = model.LinearNumExpr();
                for (int j = 0; j < data.DCSize; j++)
                {
                    if (incumbent[j] == 1)
                    {
                        consexpr.AddTerm(-1, u_j[j]);
                        constant++;
                    }
                    else
                    {
                        consexpr.AddTerm(1, u_j[j]);
                    }
                }
                model.AddLe(consexpr, radiusk - constant, "left branch");
            }
            else//right
            {
                if(incumbent.ToList().Sum() > 0)
                {
                    int constant = 0;
                    ILinearNumExpr consexpr = model.LinearNumExpr();
                    for (int j = 0; j < data.DCSize; j++)
                    {
                        if (incumbent[j] == 1)
                        {
                            consexpr.AddTerm(-1, u_j[j]);
                            constant++;
                        }
                        else
                        {
                            consexpr.AddTerm(1, u_j[j]);
                        }
                    }
                    model.AddGe(consexpr, radiusk + 1 - constant, "right branch");
                }                
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

        /// <summary>
        /// benders dual subproblem
        /// is used for the ccg master problem,
        /// which means that coefficient \xi is given
        /// and return the valid benders cut for ccg master problem
        /// </summary>
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
                delta_r[r] = model.NumVar(0, float.MaxValue, NumVarType.Float, $"delta_{r}");
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
            
            //first part            
            ILinearNumExpr subobj = model.LinearNumExpr();
            for (int r = 0; r < data.pathSize; r++)
            {
                subobj.AddTerm(-y_sol[data.s_r[r]]*(1 - u_sol[data.s_r[r]]), delta_r[r]);
                
            }
            for (int i = 0; i < data.nodeSize; i++)
            {
                subobj.AddTerm(1, v_i[i]);
            }
            if (data.Capacitystrategy)
            {
                for (int a = 0; a < data.linkSize; a++)
                {
                    pi_a[a] = model.NumVar(0, float.MaxValue, NumVarType.Float, $"pi_{a}");
                    number_of_var++;
                }

                for (int a = 0; a < data.linkSize; a++)
                {
                    subobj.AddTerm(-data.linklist[a].LinkCapacity, pi_a[a]);
                }
            }
            model.AddMaximize(subobj);

            int cons_0 = 1;
            if (cons_0 == 1)
            {                
                for (int r = 0; r < data.pathSize; r++)
                {
                    ILinearNumExpr consexpr = model.LinearNumExpr();
                    consexpr.AddTerm(-data.h_i[data.e_r[r]], w_j[data.s_r[r]]);

                    consexpr.AddTerm(1, v_i[data.e_r[r]]);

                    consexpr.AddTerm(-1, delta_r[r]);

                    if (data.Capacitystrategy)
                    {
                        for (int l = 0; l < data.linkSize; l++)
                        {
                            if (data.pathlist[r].Contains(l))
                            {
                                consexpr.AddTerm(-data.h_i[data.e_r[r]], pi_a[l]);
                            }
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
                    IRange constraint = model.AddEq(w_j[j], data.c_j[j], $"Z_Dual_{j}");
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
        /// <summary>
        /// generate upper bound of delta
        /// </summary>
        /// <param name="y_sol"></param>
        /// <param name="u_sol"></param = y_sol>
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
                delta_r[r] = model.NumVar(0, float.MaxValue, NumVarType.Float, $"delta_{r}");
                number_of_var++;
            }
            for (int i = 0; i < data.nodeSize; i++)
            {
                v_i[i] = model.NumVar(0, float.MaxValue, NumVarType.Float, $"v_{i}");
                number_of_var++;
            }
            for (int j = 0; j < data.DCSize; j++)
            {
                w_j[j] = model.NumVar(0, float.MaxValue, NumVarType.Float, $"w_{j}");
                number_of_var++;
            }

            //first part            
            ILinearNumExpr subobj = model.LinearNumExpr();
            for (int r = 0; r < data.pathSize; r++)
            {
                subobj.AddTerm(-y_sol[data.s_r[r]] * (1 - u_sol[data.s_r[r]]), delta_r[r]);
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
                    double sumofdelta = 0;
                    for (int j = 0; j < data.DCSize; j++)
                    {
                        sumofdelta += data.h_i[data.e_r[r]] * data.pathtoStartDC[r, j] * data.c_j[j];                        
                    }

                    for (int i = 0; i < data.nodeSize; i++)
                    {
                        consexpr.AddTerm(data.pathtoEndnode[r, i], v_i[i]);
                    }
                    consexpr.AddTerm(-1, delta_r[r]);

                    model.AddLe(consexpr, data.d_r[r] * data.h_i[data.e_r[r]] + sumofdelta, $"X_dual[{r}]");
                }
            }

            int cons_1 = 0;
            if (cons_1 == 1)
            {
                for (int j = 0; j < data.DCSize; j++)
                {
                    IRange constraint = model.AddLe(w_j[j], data.c_j[j], $"Z_Dual_{j}");
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
        /// <summary>
        /// generate pareto cut after identifying optimality
        /// </summary>
        /// <param name="y_solution"></param>        
        /// <param name="xI_solution"></param>
        /// <param name="objcostofdualsub"></param>
        /// <param name="y0_solution"></param>        
        /// <param name="BDM"></param>
        /// <returns></returns>
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
                    subobj.AddTerm(-y0_sol[data.s_r[r]] * (1 - u0_sol[data.s_r[r]]), delta_r[r]);
                }
                for (int i = 0; i < data.nodeSize; i++)
                {
                    subobj.AddTerm(1, v_i[i]);
                }

                original_objf.ClearExpr();
                original_objf.Expr = subobj;

            }
            // can obtain the same dual values at least
            int cons = 1;
            if (cons == 1)
            {
                ILinearNumExpr subobj = model.LinearNumExpr();
                for (int r = 0; r < data.pathSize; r++)
                {
                    subobj.AddTerm(-y_sol[data.s_r[r]] * (1 - u_sol[data.s_r[r]]), delta_r[r]);
                }
                for (int i = 0; i < data.nodeSize; i++)
                {
                    subobj.AddTerm(1, v_i[i]);
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
                    double[] v_i_values = new double[data.nodeSize];//>0
                    double[] delta_r_values = new double[data.pathSize];//>0

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

                    ILinearNumExpr extrem_point_exp = model.LinearNumExpr();

                    extrem_point_exp.AddTerm(1, BDM.omega);

                    //first part: y-related
                    for (int r = 0; r < data.pathSize; r++)
                    {
                        extrem_point_exp.AddTerm(delta_r_values[r] * (1 - u_sol[data.s_r[r]]), BDM.y_j[data.s_r[r]]);
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
        /// <summary>
        /// generate pareto cut after identifying optimality
        /// </summary>
        /// <param name="y_solution"></param>        
        /// <param name="xI_solution"></param>
        /// <param name="objcostofdualsub"></param>
        /// <param name="y0_solution"></param>        
        /// <param name="BDM"></param>
        /// <returns></returns>
        public IRange GenParetoCutforCCGMaster(int[] y_sol, int[] u_sol, double objcostofdualsub,
            int[] y0_sol, int[] u0_sol, CCGMaster BDM)
        {
            #region
            IObjective original_objf = model.GetObjective();

            int objflag = 1;
            if (objflag == 1)
            {
                ILinearNumExpr subobj = model.LinearNumExpr();
                for (int r = 0; r < data.pathSize; r++)
                {
                    subobj.AddTerm(-y0_sol[data.s_r[r]] * (1 - u0_sol[data.s_r[r]]), delta_r[r]);
                }
                for (int i = 0; i < data.nodeSize; i++)
                {
                    subobj.AddTerm(1, v_i[i]);
                }

                original_objf.ClearExpr();
                original_objf.Expr = subobj;

            }
            // can obtain the same dual values at least
            int cons = 1;
            if (cons == 1)
            {
                ILinearNumExpr subobj = model.LinearNumExpr();
                for (int r = 0; r < data.pathSize; r++)
                {
                    subobj.AddTerm(-y_sol[data.s_r[r]] * (1 - u_sol[data.s_r[r]]), delta_r[r]);
                }
                for (int i = 0; i < data.nodeSize; i++)
                {
                    subobj.AddTerm(1, v_i[i]);
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
                    double[] v_i_values = new double[data.nodeSize];//>0
                    double[] delta_r_values = new double[data.pathSize];//>0

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

                    ILinearNumExpr extrem_point_exp = model.LinearNumExpr();

                    extrem_point_exp.AddTerm(1, BDM.omega);

                    //first part: y-related
                    for (int r = 0; r < data.pathSize; r++)
                    {
                        extrem_point_exp.AddTerm(delta_r_values[r] * (1 - u_sol[data.s_r[r]]), BDM.y_j[data.s_r[r]]);
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
        /// <summary>
        /// 
        /// </summary>
        /// <param name="upd_y_solution"></param updated y>        
        /// <param name="upd_xi_j"></param updated xi>
        public void ResetBDDSubObj(int[] upd_y_sol,  int[] upd_u_j)
        {
            #region
            IObjective ccgsubobj = model.GetObjective();
            ccgsubobj.ClearExpr();
            ILinearNumExpr subobj = model.LinearNumExpr();
            for (int r = 0; r < data.pathSize; r++)
            {
                subobj.AddTerm(-upd_y_sol[data.s_r[r]] * (1 - upd_u_j[data.s_r[r]]), delta_r[r]);
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

        /// <summary>
        /// transfer the current to a new model for the consideration 
        /// of branch and check scheme when facing a new schenario
        /// </summary>
        public SCNRMaster BranchandCutForWorstCsenario(int[] y_sol)
        {
            #region            
            SCNRMaster ssmaster = new SCNRMaster(y_sol);
            xI_sols_pool = new List<int[]>();//save suboptimal scenarios
            Random rand = new Random();

            //running the final MILP model
            ssmaster.GenScenarioMaster();//generate master problem

            List<DataStructure> p_predefinedcut_list = new List<DataStructure>();

            if (data.localBranchStrategy)
            {
                List<int> opendcenterlist = new List<int>();
                List<double> correspondtransportcost = new List<double>();
                //find all possible routes with their cost according to open y
                for (int j = 0; j < data.DCSize; j++)
                {
                    if (y_sol[j] == 1)
                    {
                        
                        List<double> costlist = new List<double>();
                        for (int r = 0; r < data.pathSize; r++)
                        {
                            if(data.s_r[r] == j)
                            {
                                costlist.Add(data.d_r[r]);
                            }                            
                        }
                        opendcenterlist.Add(j);
                        correspondtransportcost.Add(costlist.Average());
                    }
                }

                //determine the worst one
                int[] worstcase = new int[data.DCSize];
                List<int> tempdcenter = opendcenterlist.ToList();
                List<double> temptransport = correspondtransportcost.ToList();
                List<int> worstlist = new List<int>();
                int DCcounter = 0;
                while (true)
                {
                    if (tempdcenter.Count == 0)
                        break;
                    if (DCcounter >= data.max_dstroyed_DCs)
                        break;

                    int ind = temptransport.IndexOf(temptransport.Min());

                    worstcase[tempdcenter[ind]] = 1;
                    worstlist.Add(tempdcenter[ind]);

                    temptransport.RemoveAt(ind); tempdcenter.RemoveAt(ind);
                    DCcounter++;
                }

                int cutcounter = 0;
                double bestobj = 0;
                SCNRSub modifiedSub = new SCNRSub(y_sol, worstcase, data);
                modifiedSub.GenScenarioSubproblem();
                modifiedSub.model.Solve();
                double benchmarkobj = modifiedSub.model.GetObjValue();
                bestobj = benchmarkobj;
                DataStructure ds = new DataStructure(data);

                for (int r = 0; r < data.pathSize; r++)
                {
                    ds.x_r_val[r] = modifiedSub.model.GetValue(modifiedSub.x_r[r]);
                }
                ds.righthandsideU = benchmarkobj;
                p_predefinedcut_list.Add(ds);

                bool termination = true;
                Stopwatch sw = new Stopwatch();
                double totallocalbranchseconds = 0;

                while (termination)
                {
                    sw.Start();
                    int[] tempworstcase = worstcase.ToArray();

                    if (cutcounter < 5)
                    {
                        tempworstcase[worstlist.Last()] = 0;
                        int localoperator = rand.Next(0, tempdcenter.Count);
                        tempworstcase[tempdcenter[localoperator]] = 1;
                    }
                    else if (cutcounter < 10)
                    {
                        List<int> p_tempcase = worstlist.ToList();

                        for (int i = 0; i < 2; i++)
                        {
                            int index = rand.Next(0, p_tempcase.Count);
                            tempworstcase[p_tempcase[index]] = 0;
                            p_tempcase.RemoveAt(index);

                            int localoperator = rand.Next(0, tempdcenter.Count);
                            tempworstcase[tempdcenter[localoperator]] = 1;

                            if (p_tempcase.Count == 0)
                                break;
                        }

                    }
                    else
                    {
                        break;
                    }

                    modifiedSub.ResetOBJaRHS_XI(y_sol, tempworstcase, data.big_M_r);
                    modifiedSub.model.Solve();
                    double currentobj = modifiedSub.model.GetObjValue();

                    if (currentobj > bestobj)
                    {
                        bestobj = currentobj;
                        ds = new DataStructure(data);

                        for (int r = 0; r < data.pathSize; r++)
                        {
                            ds.x_r_val[r] = modifiedSub.model.GetValue(modifiedSub.x_r[r]);
                        }
                        ds.righthandsideU = currentobj;
                        p_predefinedcut_list.Add(ds);
                        cutcounter = 0;
                    }
                    else
                    {
                        cutcounter++;
                    }
                    sw.Stop();
                    totallocalbranchseconds += sw.ElapsedMilliseconds / 1000;
                    if (totallocalbranchseconds > 3)
                        break;
                }
            }

            ssmaster.Addpredefinedcuts(p_predefinedcut_list);

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
        public void generateLshapedbound()
        {
            #region
            CCGMaster ccgmst = new CCGMaster();
            ccgmst.InitializeCCGMaster();
            ccgmst.u_l_j.Add(new int[data.DCSize]);
            ccgmst.GenCCGMasterproblem(1);
            ccgmst.model.Solve();

            data.LShapedbound = ccgmst.model.GetValue(ccgmst.omega);
           
            ccgmst = null;

            #endregion
        }
        public void generateParetoPoint()
        {
            #region
            CCGMaster ccgmst = new CCGMaster();
            //generate input node for pareto cut
            ccgmst.InitializeCCGMaster();
            ccgmst.GenParetoCorePoint();
            ccgmst.model.Solve();

            data.pareto_y_sol = new int[data.DCSize];
            
            for (int j = 0; j < data.DCSize; j++)
            {
                double y_value = ccgmst.model.GetValue(ccgmst.y_j[j]);
                if (y_value > 0.5)
                {
                    data.pareto_y_sol[j] = 1;
                }                
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
        /// <summary>
        /// transfer the current to a new model for the consideration 
        /// of branch and check scheme when facing a new schenario
        /// </summary>
        public BendersMaster BranchandCut()
        {
            #region

            Program.SolutionIteration.WriteLine("findworsttime, subproblemtime,multicutstime,paretotime,worst_scenario_cost,best objective cost");
            Program.g_comparativeIndicators.WriteLine("# feasol, # optsol, best integer, best bound, Gap, CPU time");
            BendersMaster bdmaster = new BendersMaster();
            CCGMaster CCGMasterProblem = new CCGMaster();
            //generate lshapedbound
            generateLshapedbound();

            Stopwatch optimize_BAC = new Stopwatch();
            CCGMasterProblem.YCoeff_list = new List<DataStructure>();

            if (data.stablilization_at_rootnode)
            {
                //using stabilized node to accelerate solution process
                int total_sets = 3;
                double[] alpha_list = new double[total_sets];
                double[] lambda_list = new double[total_sets];
                //generate multiple parameters for parallel computing on stabilization
                for (int m = 0; m < total_sets; m++)
                {
                    double lambda_val = data.Para_STB_Lambada + 0.1 * m;
                    double alpha_val = data.Para_STB_Alpha - 0.1 * m;

                    alpha_list[m] = alpha_val;
                    lambda_list[m] = lambda_val;
                }

                object lockObject = new object();//protect the file
                var options = new ParallelOptions() { MaxDegreeOfParallelism = Environment.ProcessorCount };


                Parallel.For(0, total_sets, options, (i, state) =>
                {
                    CCGMaster ccgmp = new CCGMaster();

                    ccgmp.BendersDec(lambda_list[i], data.Para_STB_Epsilon, alpha_list[i]);

                    lock (lockObject)
                    {
                        for (int f = 0; f < ccgmp.Temp_YCoeff_list.Count; f++)
                        {
                            CCGMasterProblem.YCoeff_list.Add(ccgmp.Temp_YCoeff_list[f]);
                        }
                    }
                });
            }

            generateParetoPoint();//generate pareto point

            //generate benders master problem
            bdmaster.GenBDMasterproblem();

            //Add cuts obtained from stabilized step
            bdmaster.AddNewValidCuts(CCGMasterProblem.YCoeff_list);
            CCGMasterProblem = null;
            // set time limitaion
            bdmaster.model.SetParam(Cplex.Param.TimeLimit, data.TL);
            // Set up the cut callback to be used for separating Benders' cuts
            bdmaster.model.SetParam(Cplex.Param.Preprocessing.Presolve, false);
            bdmaster.model.SetParam(Cplex.Param.Threads, 1);
            // Turn on traditional search for use with control callbacks
            bdmaster.model.SetParam(Cplex.Param.MIP.Strategy.Search, Cplex.MIPSearch.Traditional);
            bdmaster.model.Use(new BendersLazyConsCallback(bdmaster, data));
            //bdmaster.model.Use(new BendersHeuristicCallback(bdmaster, data));

            bdmaster.model.SetOut(Program.TWoutput);
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

            Cplex.Status solvests = bdmaster.model.GetStatus();

            Program.SolutionIteration.WriteLine($"Total B&C time,{totalCPUtime}");
            Program.SolutionIteration.WriteLine($"Optimal obj,{optobj}");
            Program.SolutionIteration.WriteLine($"Best Linear bound,{bestbound}");

            return bdmaster;

            #endregion
        }
        /// <summary>
        /// transfer the current to a new model for the consideration 
        /// of branch and check scheme when facing a new schenario
        /// </summary>
        public CCGMaster BranchandCheckforCCGIteration(CCGMaster ccgmaster, List<int[]> selectedScenarios)
        {
            #region

            Program.SolutionIteration.WriteLine("findworsttime, subproblemtime,multicutstime,paretotime,worst_scenario_cost,best objective cost");
            Program.g_comparativeIndicators.WriteLine("# feasol, # optsol, best integer, best bound, Gap, CPU time");

            CCGMaster CCGMasterProblem = new CCGMaster();
            //generate lshapedbound
            generateLshapedbound();

            Stopwatch optimize_BAC = new Stopwatch();
            CCGMasterProblem.YCoeff_list = new List<DataStructure>();

            if (data.stablilization_at_rootnode)
            {
                //using stabilized node to accelerate solution process
                int total_sets = 3;
                double[] alpha_list = new double[total_sets];
                double[] lambda_list = new double[total_sets];
                //generate multiple parameters for parallel computing on stabilization
                for (int m = 0; m < total_sets; m++)
                {
                    double lambda_val = data.Para_STB_Lambada + 0.1 * m;
                    double alpha_val = data.Para_STB_Alpha - 0.1 * m;

                    alpha_list[m] = alpha_val;
                    lambda_list[m] = lambda_val;
                }

                object lockObject = new object();//protect the file
                var options = new ParallelOptions() { MaxDegreeOfParallelism = Environment.ProcessorCount };


                Parallel.For(0, total_sets, options, (i, state) =>
                {
                    CCGMaster ccgmp = new CCGMaster();

                    ccgmp.BendersDec(lambda_list[i], data.Para_STB_Epsilon, alpha_list[i]);

                    lock (lockObject)
                    {
                        for (int f = 0; f < ccgmp.Temp_YCoeff_list.Count; f++)
                        {
                            CCGMasterProblem.YCoeff_list.Add(ccgmp.Temp_YCoeff_list[f]);
                        }
                    }
                });
            }

            generateParetoPoint();//generate pareto point

            //Add cuts obtained from stabilized step
            ccgmaster.AddNewValidCuts(CCGMasterProblem.YCoeff_list);
            CCGMasterProblem = null;
            // set time limitaion
            ccgmaster.model.SetParam(Cplex.Param.TimeLimit, data.TL);
            // Set up the cut callback to be used for separating Benders' cuts
            ccgmaster.model.SetParam(Cplex.Param.Preprocessing.Presolve, false);
            ccgmaster.model.SetParam(Cplex.Param.Threads, 1);
            // Turn on traditional search for use with control callbacks
            ccgmaster.model.SetParam(Cplex.Param.MIP.Strategy.Search, Cplex.MIPSearch.Traditional);
            ccgmaster.model.Use(new BLCCforCCGMaster(ccgmaster, data, selectedScenarios));
            ccgmaster.model.Use(new BendersUserCutCallback(ccgmaster, data));

            ccgmaster.model.SetOut(Program.TWoutput);
            optimize_BAC.Start();
            ccgmaster.model.Solve();
            optimize_BAC.Stop();

            double optobj = ccgmaster.model.GetObjValue();
            double bestbound = ccgmaster.model.GetBestObjValue();
            double totalCPUtime = optimize_BAC.ElapsedMilliseconds / 1000;
            double relativeGap = ccgmaster.model.GetMIPRelativeGap();
            int totalSolutionCount = ccgmaster.model.GetSolnPoolNsolns();
            int Nofeasols = 0;
            int Nooptsols = 0;

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
                        if (ccgmaster.model.GetValue(ccgmaster.y_j[j], n) > 0.5)
                        {
                            y_val[j] = 1;
                        }
                        omega_val = ccgmaster.model.GetValue(ccgmaster.omega, n);
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

            Cplex.Status solvests = ccgmaster.model.GetStatus();

            Program.SolutionIteration.WriteLine($"Total B&C time,{totalCPUtime}");
            Program.SolutionIteration.WriteLine($"Optimal obj,{optobj}");
            Program.SolutionIteration.WriteLine($"Best Linear bound,{bestbound}");

            return ccgmaster;

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

        public void ParallelCCG()
        {
            #region solution process
            //generate master problem and subproblem
            Program.g_CCGiteration.WriteLine("IterNumber, UB, LB, etaval, suproblemtime, mastertime");
            Program.g_CPLEXResults.WriteLine("# feasol, # optsol, best integer, best bound, Gap, CPU time");

            List<Path> feedbacktomaster = new List<Path>();
            double eta_value = 0, mastercost = 0, UB = float.MaxValue, LB = float.MinValue;
            bool isparrallel = true;
            int cutting_flag = 0;// add the cutting plane
            int[] y_solution = new int[data.DCSize];

            Stopwatch subpro = new Stopwatch();
            Stopwatch masterpro = new Stopwatch();
            Stopwatch totalProcedure = new Stopwatch();
            double mastertime = 0;
            double subtime = 0;

            double mastertime_iter = 0;
            double subtime_iter = 0;

            totalProcedure.Start();
            double availabletime = data.TL;
            masterproblem.InitializeCCGMaster();
            masterproblem.model.SetParam(Cplex.Param.Threads, 1);
            masterproblem.model.SetParam(Cplex.Param.TimeLimit, availabletime);
            int number_of_opens = 0;
            for (int j = 0; j < data.DCSize; j++)
            {
                y_solution[j] = 1;
                number_of_opens++;
            }

            bool feasible_flag = false;
            double worst_scenario_cost = float.MinValue;

            subpro = new Stopwatch();
            subpro.Start();

            int[] worst_u_j = new int[data.DCSize];
            
            List<int[]> scenarios = new List<int[]>();
            List<Cplex> subproblemlist = new List<Cplex>();
            List<CCGSub> ccgsublist = new List<CCGSub>();

            int worstindex = 0;//for outputing results
            int number_of_scenarios = 0;
            // Find the worst scenario in the context of the current solution, i.e., two cases
            //case 1: the number of open DCents less than or equal to the predefined threshold
            worst_scenario_cost = float.MinValue;
            object lockObject = new object();//protect the file
            
            feasible_flag = true;//indicate whether the subproblem is feasible or not
            List<double> objlist = new List<double>();
            List<int[]> cuttingflaglist = new List<int[]>();

            List<int> open_DCent_No_list = new List<int>();
            for (int j = 0; j < data.DCSize; j++)
            {
                if (y_solution[j] == 1)
                {
                    open_DCent_No_list.Add(j);
                }
            }

            // All scenarios
            List<List<int>> all_scenarios = new List<List<int>>();

            //no decents is destroyed                
            //enumarate all scenarios
            List<List<int>> temp_scenarios = Combination.Combine(open_DCent_No_list, data.max_dstroyed_DCs);

            for (int l = 0; l < temp_scenarios.Count; l++)
            {
                List<int> tem_destroy_solution = temp_scenarios[l].ToList();

                all_scenarios.Add(tem_destroy_solution.ToList());
            }

            number_of_scenarios = all_scenarios.Count;

            for (int l = 0; l < all_scenarios.Count; l++)
            {
                List<int> temp_solution = all_scenarios[l].ToList();

                int[] Dcent_state = new int[data.DCSize];

                for (int j = 0; j < data.DCSize; j++)
                {
                    Dcent_state[j] = 0;//initialization;                        
                }

                for (int j = 0; j < temp_solution.Count; j++)
                {
                    int dcentno = temp_solution[j];
                    Dcent_state[dcentno] = 1;//diruption                       
                }
                subproblem = new CCGSub();
                subproblem.GenCCGSubproblem(y_solution, Dcent_state);
                //subproblem.model.ExportModel("CCGsubproblem.lp");
                //collect subproblem                  
                subproblemlist.Add(subproblem.model);
                ccgsublist.Add(subproblem);
                //collect scenarios
                scenarios.Add(Dcent_state);
                if (!isparrallel)
                {
                    if (subproblem.model.Solve())
                    {
                        if (subproblem.model.GetObjValue() > worst_scenario_cost)
                        {
                            worst_scenario_cost = subproblem.model.GetObjValue();
                            worst_u_j = Dcent_state.ToArray();
                        }
                    }
                    else
                    {
                        feasible_flag = false;
                        worst_u_j = Dcent_state.ToArray();
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
                        lock (lockObject)
                        {
                            infeasibleFound = true;
                            worst_u_j = scenarios[i].ToArray();
                        }
                        state.Stop();
                    }

                });
                if (!infeasibleFound)
                {
                    for (int s = 0; s < subproblemlist.Count; s++)
                    {
                        double objcost = subproblemlist[s].GetObjValue();
                        if (objcost > worst_scenario_cost)
                        {
                            worst_scenario_cost = objcost;
                            worst_u_j = scenarios[s].ToArray();
                            worstindex = s;
                        }

                        if (data.MultiplescenariosforCCG)
                        {
                            //data.maxMultiCuts + 1 includes the worst scenario
                            if (objlist.Count < data.maxMultiCuts + 1)
                            {
                                objlist.Add(worst_scenario_cost);
                                cuttingflaglist.Add(worst_u_j);
                            }
                            else
                            {
                                double minobj = objlist.Min();
                                int indx = objlist.IndexOf(minobj);

                                if (minobj < objcost)
                                {
                                    objlist.RemoveAt(indx);
                                    cuttingflaglist.RemoveAt(indx);

                                    objlist.Add(objcost);
                                    cuttingflaglist.Add(scenarios[s].ToArray());
                                }
                            }
                        }
                        else
                        {
                            cuttingflaglist = new List<int[]>();
                            cuttingflaglist.Add(worst_u_j);
                        }
                        //subproblemlist[s].End();
                    }
                }
                //if (!infeasibleFound)
                //{
                //    for (int s = 0; s < subproblemlist.Count; s++)
                //    {
                //        double objcost = subproblemlist[s].GetObjValue();                        
                //        if (objcost > worst_scenario_cost)
                //        {
                //            worst_scenario_cost = objcost;
                //            worst_u_j = scenarios[s].ToArray();
                //            worstindex = s;
                //        }
                //        //subproblemlist[s].End();
                //    }
                //}
            }
            subpro.Stop();

            subtime_iter = subpro.ElapsedMilliseconds;
            subtime += subtime_iter;

            Program.g_CCGiteration.WriteLine($"{iter_l},{UB}, {LB},{eta_value}, {subtime_iter}, {mastertime_iter}");
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
            
            //cuttingflaglist.Add(worst_u_j); objlist.Add(worst_scenario_cost);

            double eps = 0.001;
            iter_l = -1;
            double tickcountstart = Environment.TickCount;
            double totalrunningtime = 0;

            while (UB - LB > eps && totalrunningtime <= data.TL)
            {
                iter_l++;

                masterpro = new Stopwatch();
                masterpro.Start();

                if (cuttingflaglist.Count != 0)// existing a worst scenario
                {
                    cutting_flag = 1;//add the cutting plane
                    for (int s = 0; s < cuttingflaglist.Count; s++)
                    {
                        int[] candidateuj = cuttingflaglist[s];

                        masterproblem.u_l_j.Add(candidateuj);
                        masterproblem.GenCCGMasterproblem(cutting_flag);
                        masterproblem.model.SetParam(Cplex.Param.TimeLimit, availabletime);
                    }
                }
                else//now no feasible solutions
                {
                    cutting_flag = 0;//add the cutting plane

                    masterproblem.u_l_j.Add(worst_u_j);
                    masterproblem.GenCCGMasterproblem(cutting_flag);
                    masterproblem.model.SetParam(Cplex.Param.TimeLimit, availabletime);
                }

                masterproblem.model.Solve();
                //masterproblem.model.ExportModel("masterproblemmodel.lp");
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
                
                cuttingflaglist = new List<int[]>();
                objlist = new List<double>();

                subpro = new Stopwatch();
                subpro.Start();

                worst_u_j = new int[data.DCSize];
                scenarios = new List<int[]>();

                // Find the worst scenario in the context of the current solution, i.e., two cases
                //case 1: the number of open DCents less than or equal to the predefined threshold
                worst_scenario_cost = float.MinValue;
                feasible_flag = true;//indicate whether the subproblem is feasible or not
                subproblemlist = new List<Cplex>();
                ccgsublist = new List<CCGSub>();

                open_DCent_No_list = new List<int>();
                for (int j = 0; j < data.DCSize; j++)
                {
                    if (y_solution[j] == 1)
                    {
                        open_DCent_No_list.Add(j);
                    }
                }
                masterpro.Stop();

                mastertime_iter = masterpro.ElapsedMilliseconds;
                mastertime += mastertime_iter;

                // All scenarios
                all_scenarios = new List<List<int>>();

                //no decents is destroyed                
                //enumarate all scenarios
                temp_scenarios = Combination.Combine(open_DCent_No_list, data.max_dstroyed_DCs);

                for (int l = 0; l < temp_scenarios.Count; l++)
                {
                    List<int> tem_destroy_solution = temp_scenarios[l].ToList();
                    all_scenarios.Add(tem_destroy_solution.ToList());
                }
                number_of_scenarios = all_scenarios.Count;

                for (int l = 0; l < all_scenarios.Count; l++)
                {
                    List<int> temp_solution = all_scenarios[l].ToList();

                    int[] Dcent_state = new int[data.DCSize];

                    for (int j = 0; j < data.DCSize; j++)
                    {
                        Dcent_state[j] = 0;//initialization;                        
                    }

                    for (int j = 0; j < temp_solution.Count; j++)
                    {
                        int dcentno = temp_solution[j];
                        Dcent_state[dcentno] = 1;//diruption                       
                    }

                    subproblem = new CCGSub();
                    subproblem.GenCCGSubproblem(y_solution, Dcent_state);
                    subproblem.model.SetParam(Cplex.Param.Threads, 1);
                    //subproblem.model.ExportModel("CCGsubproblem.lp");
                    //collect subproblem                  
                    subproblemlist.Add(subproblem.model);
                    ccgsublist.Add(subproblem);
                    //collect scenarios
                    scenarios.Add(Dcent_state);
                    if (!isparrallel)
                    {
                        if (subproblem.model.Solve())
                        {
                            if (subproblem.model.GetObjValue() > worst_scenario_cost)
                            {
                                worst_scenario_cost = subproblem.model.GetObjValue();
                                worst_u_j = Dcent_state.ToArray();

                            }
                        }
                        else
                        {
                            feasible_flag = false;
                            worst_u_j = Dcent_state.ToArray();
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
                            lock (lockObject)
                            {
                                infeasibleFound = true;
                                feasible_flag = false;
                                worst_u_j = scenarios[i].ToArray();
                                worstindex = i;
                            }
                            state.Stop();
                        }

                    });
                    if (!infeasibleFound)
                    {
                        for (int s = 0; s < subproblemlist.Count; s++)
                        {
                            double objcost = subproblemlist[s].GetObjValue();
                            if (objcost > worst_scenario_cost)
                            {
                                worst_scenario_cost = objcost;
                                worst_u_j = scenarios[s].ToArray();
                                worstindex = s;                                
                            }

                            if (data.MultiplescenariosforCCG)
                            {
                                //data.maxMultiCuts + 1 includes the worst scenario
                                if (objlist.Count < data.maxMultiCuts + 1)
                                {
                                    objlist.Add(worst_scenario_cost);
                                    cuttingflaglist.Add(worst_u_j);
                                }
                                else
                                {
                                    double minobj = objlist.Min();
                                    int indx = objlist.IndexOf(minobj);

                                    if (minobj < objcost)
                                    {
                                        objlist.RemoveAt(indx);
                                        cuttingflaglist.RemoveAt(indx);

                                        objlist.Add(objcost);
                                        cuttingflaglist.Add(scenarios[s].ToArray());
                                    }
                                }
                            }
                            else
                            {
                                cuttingflaglist = new List<int[]>();
                                cuttingflaglist.Add(worst_u_j);
                            }
                            //subproblemlist[s].End();
                        }
                    }                    
                }
                subpro.Stop();

                subtime_iter = subpro.ElapsedMilliseconds;
                subtime += subtime_iter;

                //find the worst scenario
                if (feasible_flag)
                {
                    //find the worst scenario
                    UB = Math.Min(UB, worst_scenario_cost + mastercost - eta_value);
                }

                Program.g_CCGiteration.WriteLine($"{iter_l},{UB}, {LB},{eta_value}, {subtime_iter}, {mastertime_iter}");
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
                totalrunningtime = (Environment.TickCount - tickcountstart) / 1000;
                availabletime = data.TL - totalrunningtime;
            }

            totalProcedure.Stop();

            double optobj = masterproblem.model.GetObjValue();
            double bestbound = masterproblem.model.GetBestObjValue();
            double totalCPUtime = totalProcedure.ElapsedMilliseconds / 1000;
            double relativeGap = masterproblem.model.GetMIPRelativeGap();

            if (UB - LB > eps && totalrunningtime > data.TL)
            {
                optobj = UB; bestbound = LB; relativeGap = (UB - LB) / UB * 100;
            }

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

            Program.g_CCGiteration.WriteLine();
            Program.g_CCGiteration.WriteLine($"Totalmastertime,{mastertime}");
            Program.g_CCGiteration.WriteLine($"Totalsubproblemtime,{subtime}");

            solution.write_solution(masterproblem);            
            solution.write_subproblemsolution(ccgsublist[worstindex]);

            Console.WriteLine("UB:{0}, LB:{1}", UB, LB);
            #endregion
        }
        public void CCGwBAC()
        {
            #region solution process
            //generate master problem and subproblem
            Program.g_CCGiteration.WriteLine("IterNumber, UB, LB, etaval, suproblemtime, mastertime");
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
            int cutting_flag = 0;// add the cutting plane
            int[] y_solution = new int[data.DCSize];
            subproblem = new CCGSub();
            subproblem.GenCCGSubproblem(y_solution, new int[data.DCSize]);
            masterproblem.InitializeCCGMaster();
            masterproblem.model.SetParam(Cplex.Param.Threads, 1);
            masterproblem.model.SetParam(Cplex.Param.TimeLimit, availabletime);

            int number_of_opens = 0;
            for (int j = 0; j < data.DCSize; j++)
            {

                y_solution[j] = 1;
                number_of_opens++;
            }            

            bool feasible_flag = false;
            double worst_scenario_cost = float.MinValue;

            subpro = new Stopwatch();
            subpro.Start();
            double tickcountstart = Environment.TickCount;
            double totalrunningtime = 0;
            int[] worst_u_j = new int[data.DCSize];

            CACGBD colBDGen = new CACGBD();

            SCNRMaster submaster = colBDGen.BranchandCutForWorstCsenario(y_solution);

            Cplex.Status cur_status = submaster.model.GetStatus();

            if (!submaster.feasiblestatus)
            {
                worst_u_j = submaster.worstScenarioSolution.ToArray();
            }
            else
            {
                if (cur_status == Cplex.Status.Optimal || cur_status == Cplex.Status.Feasible)
                    feasible_flag = true;
                else if (cur_status == Cplex.Status.Infeasible)
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
            }
            
            subpro.Stop();

            subtime_iter = subpro.ElapsedMilliseconds;
            subtime += subtime_iter;

            Program.g_CCGiteration.WriteLine($"{iter_l},{UB}, {LB},{eta_value}, {subtime_iter}, {mastertime_iter}");
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
            cuttingflaglist.Add(worst_u_j);

            double eps = 0.0001;
            iter_l = -1;

            while (UB - LB > eps * UB && totalrunningtime <= data.TL)
            {
                iter_l++;
                masterpro = new Stopwatch();
                masterpro.Start();

                if (cuttingflaglist.Count != 0)// existing a worst scenario
                {
                    cutting_flag = 1;//add the cutting plane
                    for (int s = 0; s < cuttingflaglist.Count; s++)
                    {
                        int[] candidateuj = cuttingflaglist[s];

                        masterproblem.u_l_j.Add(candidateuj);
                        masterproblem.GenCCGMasterproblem(cutting_flag);
                        masterproblem.model.SetParam(Cplex.Param.TimeLimit, availabletime);
                    }
                }
                else//now no feasible solutions
                {
                    cutting_flag = 0;//add the cutting plane

                    masterproblem.u_l_j.Add(worst_u_j);
                    masterproblem.GenCCGMasterproblem(cutting_flag);
                    masterproblem.model.SetParam(Cplex.Param.TimeLimit, availabletime);
                }
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

                masterpro.Stop();

                mastertime_iter = masterpro.ElapsedMilliseconds;
                mastertime += mastertime_iter;

                feasible_flag = false;
                cuttingflaglist = new List<int[]>();

                subpro = new Stopwatch();
                subpro.Start();
                worst_u_j = new int[data.DCSize];

                // Find the worst scenario in the context of the current solution
                submaster = colBDGen.BranchandCutForWorstCsenario(y_solution);

                if (!submaster.feasiblestatus)
                {
                    worst_u_j = submaster.worstScenarioSolution.ToArray();
                }
                else
                {
                    cur_status = submaster.model.GetStatus();
                    if (cur_status == Cplex.Status.Optimal || cur_status == Cplex.Status.Feasible)
                        feasible_flag = true;
                    else if (cur_status == Cplex.Status.Infeasible)
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
                    //
                    subproblem.ResetOBJaRHS_XI(y_solution, worst_u_j);
                    subproblem.model.SetParam(Cplex.Param.Preprocessing.Presolve, false);
                    subproblem.model.SetParam(Cplex.Param.Threads, 1);
                    //ccgsub.model.SetParam(Cplex.Param.RootAlgorithm, Cplex.Algorithm.Dual);
                    subproblem.model.Solve();
                    cur_status = subproblem.model.GetStatus();

                    if (cur_status == Cplex.Status.Optimal || cur_status == Cplex.Status.Feasible)
                    {
                        feasible_flag = true;
                        worst_scenario_cost = subproblem.model.GetObjValue();
                        //find the worst scenario                    
                        UB = Math.Min(UB, worst_scenario_cost + mastercost - eta_value);
                        cuttingflaglist.Add(worst_u_j);

                        if (data.MultiplescenariosforCCG)
                        {
                            for (int s = 0; s < colBDGen.xI_sols_pool.ToList().Count; s++)
                            {
                                cuttingflaglist.Add(colBDGen.xI_sols_pool.ToList()[s]);
                            }
                        }
                    }
                }
                
                subpro.Stop();

                subtime_iter = subpro.ElapsedMilliseconds;
                subtime += subtime_iter;

                Program.g_CCGiteration.WriteLine($"{iter_l},{UB}, {LB},{eta_value}, {subtime_iter}, {mastertime_iter}");
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

            double optobj = masterproblem.model.GetObjValue();
            double bestbound = masterproblem.model.GetBestObjValue();
            double totalCPUtime = totalProcedure.ElapsedMilliseconds / 1000;
            double relativeGap = masterproblem.model.GetMIPRelativeGap();

            if (UB - LB > eps && totalrunningtime > data.TL)
            {
                optobj = UB; bestbound = LB; relativeGap = (UB - LB) / UB * 100;
            }

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

            Program.g_CCGiteration.WriteLine();
            Program.g_CCGiteration.WriteLine($"Totalmastertime,{mastertime}");
            Program.g_CCGiteration.WriteLine($"Totalsubproblemtime,{subtime}");

            solution.write_solution(masterproblem);
            solution.write_subproblemsolution(subproblem);
            Console.WriteLine("UB:{0}, LB:{1}", UB, LB);
            #endregion
        }
        /// <summary>
        /// find robust plan by solving MILP
        /// </summary>
        public void CCGMILP()
        {
            #region solution process
            //generate master problem and subproblem

            Program.g_CCGiteration.WriteLine("IterNumber, UB, LB, etaval, suproblemtime, mastertime");
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

            int[] all_disrupt_u_sol = new int[data.DCSize];
            for (int j = 0; j < data.DCSize; j++)
            {
                all_disrupt_u_sol[j] = 1;
            }
            double availabletime = data.TL;
            int cutting_flag = 0;// add the cutting plane
            int[] y_solution = new int[data.DCSize];
            //initialization
            subproblem = new CCGSub();
            subproblem.GenCCGSubproblem(y_solution, new int[data.DCSize]);
            masterproblem.InitializeCCGMaster();
            masterproblem.model.SetParam(Cplex.Param.Threads, 1);
            masterproblem.model.SetParam(Cplex.Param.TimeLimit, availabletime);

            int number_of_opens = 0;
            for (int j = 0; j < data.DCSize; j++)
            {                
                y_solution[j] = 1;
                number_of_opens++;
            }

            bool feasible_flag = false;
            double worst_scenario_cost = float.MinValue;

            subpro = new Stopwatch();
            subpro.Start();
            double tickcountstart = Environment.TickCount;
            double totalrunningtime = 0;

            int[] worst_u_j = new int[data.DCSize];

            // Find the worst scenario in the context of the current solution                    
            BendersDualSub bddsub = new BendersDualSub();
            bddsub.GendeltaUBsubDual(y_solution, y_solution);
            bddsub.model.SetParam(Cplex.Param.Preprocessing.Presolve, false);
            bddsub.model.Solve();

            double[] deltaval = new double[data.pathSize];
            deltaval = bddsub.model.GetValues(bddsub.delta_r);

            dualsubproblem = new DualCCGSub();
            dualsubproblem.GenDualofCCGsub(y_solution, deltaval);

            dualsubproblem.model.SetParam(Cplex.Param.MIP.Pool.Capacity, data.maxMultiCuts); // 解池容量
            dualsubproblem.model.SetParam(Cplex.Param.MIP.Pool.Replace, 2);
            dualsubproblem.model.SetParam(Cplex.Param.MIP.Pool.Intensity, 0);
            
            dualsubproblem.model.Solve();

            Cplex.Status cur_status = dualsubproblem.model.GetStatus();
            if (cur_status == Cplex.Status.Optimal || cur_status == Cplex.Status.Feasible)
                feasible_flag = true;
            else if(cur_status == Cplex.Status.Infeasible)
            {
                Console.WriteLine("Please check!!! No feasible solutions available!!!!");
                totalrunningtime = data.TL + 1;
            }

            for (int l = 0; l < data.DCSize; l++)
            {
                var u_val = dualsubproblem.model.GetValue(dualsubproblem.u_j[l]);

                if (u_val > 0.5)
                {
                    worst_u_j[l] = 1;
                }
            }
            subpro.Stop();

            subtime_iter = subpro.ElapsedMilliseconds;
            subtime += subtime_iter;

            Program.g_CCGiteration.WriteLine($"{iter_l},{UB}, {LB},{eta_value}, {subtime_iter}, {mastertime_iter}");
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
            subproblemscenarios.Add(worst_u_j);

            double eps = 0.0001;
            iter_l = -1;

            while (UB - LB > eps * UB && totalrunningtime <= data.TL)
            {
                iter_l++;
                masterpro = new Stopwatch();
                masterpro.Start();

                if (cuttingflaglist.Count != 0)// existing a worst scenario
                {
                    cutting_flag = 1;//add the cutting plane
                    for (int s = 0; s < cuttingflaglist.Count; s++)
                    {
                        int[] candidateuj = cuttingflaglist[s];

                        masterproblem.u_l_j.Add(candidateuj);
                        masterproblem.GenCCGMasterproblem(cutting_flag);
                        masterproblem.model.SetParam(Cplex.Param.TimeLimit, availabletime);
                    }

                }
                else//now no feasible solutions
                {
                    cutting_flag = 0;//add the cutting plane

                    masterproblem.u_l_j.Add(worst_u_j);
                    masterproblem.GenCCGMasterproblem(cutting_flag);
                    masterproblem.model.SetParam(Cplex.Param.TimeLimit, availabletime);
                }

                if (data.solveCCGmasterIter)
                {
                    CACGBD solvingmaster = new CACGBD();
                    CCGMaster ccgmasteriter = solvingmaster.BranchandCheckforCCGIteration(masterproblem, subproblemscenarios);
                }
                else
                {
                    masterproblem.model.Solve();
                }                    

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
                masterpro.Stop();

                mastertime_iter = masterpro.ElapsedMilliseconds;
                mastertime += mastertime_iter;

                subpro = new Stopwatch();
                subpro.Start();

                feasible_flag = false;
                cuttingflaglist = new List<int[]>();

                worst_u_j = new int[data.DCSize];

                // Find the worst scenario in the context of the current solution                    
                bddsub.ResetBDDSubObj(y_solution, y_solution);
                bddsub.model.Solve();
                //get the upper bound of delta
                deltaval = bddsub.model.GetValues(bddsub.delta_r);

                //update deltavalue in dual of ccg subproblem
                dualsubproblem.ResetdualccgSubObjcons(y_solution, deltaval);
                dualsubproblem.model.Solve();

                Cplex.Status solvingstatus = dualsubproblem.model.GetStatus();
                for (int l = 0; l < data.DCSize; l++)
                {
                    var u_val = dualsubproblem.model.GetValue(dualsubproblem.u_j[l]);

                    if (u_val > 0.5)
                    {
                        worst_u_j[l] = 1;
                    }
                }

                if (solvingstatus == Cplex.Status.Optimal || solvingstatus == Cplex.Status.Feasible)
                {
                    //find the objective cost with worst scenario
                    subproblem.ResetOBJaRHS_XI(y_solution, worst_u_j);
                    subproblem.model.SetParam(Cplex.Param.Threads, 1);
                    subproblem.model.SetParam(Cplex.Param.Preprocessing.Presolve, false);
                    //ccgsub.model.SetParam(Cplex.Param.RootAlgorithm, Cplex.Algorithm.Dual);
                    subproblem.model.Solve();
                    cur_status = subproblem.model.GetStatus();

                    if (cur_status == Cplex.Status.Optimal || cur_status == Cplex.Status.Feasible)
                    {
                        feasible_flag = true;
                        worst_scenario_cost = subproblem.model.GetObjValue();

                        //find the worst scenario                    
                        UB = Math.Min(UB, worst_scenario_cost + mastercost - eta_value);
                        cuttingflaglist.Add(worst_u_j);

                        if (data.solveCCGmasterIter)
                        {
                            if (subproblemscenarios.Count < data.maxMultiCuts)
                            {
                                subproblemscenarios.Add(worst_u_j);
                            }
                            else
                            {
                                //First case: 
                                subproblemscenarios.RemoveAt(0);
                                subproblemscenarios.Add(worst_u_j);
                            }
                        }
                        
                        if (data.MultiplescenariosforCCG)
                        {
                            int numSols = dualsubproblem.model.GetSolnPoolNsolns();

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
                                        if (dualsubproblem.model.GetValue(dualsubproblem.u_j[j], n) > 0.5)
                                        {
                                            feasiblesolution[j] = 1;
                                        }
                                        if (feasiblesolution[j] != worst_u_j[j])
                                            opttest = false;
                                    }
                                    double objval = dualsubproblem.model.GetObjValue();

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
                                    cuttingflaglist.Add(candidatesolset[ind]);

                                    candidatesolset.RemoveAt(ind); costset.RemoveAt(ind);
                                    counter++;
                                }
                            }
                        }
                    }
                }
                else
                {
                    Console.WriteLine();
                }
                
                subpro.Stop();

                subtime_iter = subpro.ElapsedMilliseconds;
                subtime += subtime_iter;
                
                Program.g_CCGiteration.WriteLine($"{iter_l},{UB}, {LB},{eta_value}, {subtime_iter}, {mastertime_iter}");
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

            double optobj = masterproblem.model.GetObjValue();
            double bestbound = masterproblem.model.GetBestObjValue();
            double totalCPUtime = totalProcedure.ElapsedMilliseconds / 1000;
            double relativeGap = masterproblem.model.GetMIPRelativeGap();

            if(UB - LB > eps && totalrunningtime > data.TL)
            {
                optobj = UB; bestbound = LB; relativeGap = (UB - LB) / UB * 100;
            }

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

            Program.g_CCGiteration.WriteLine();
            Program.g_CCGiteration.WriteLine($"Totalmastertime,{mastertime}");
            Program.g_CCGiteration.WriteLine($"Totalsubproblemtime,{subtime}");

            solution.write_solution(masterproblem);
            solution.write_subproblemsolution(subproblem);
            Console.WriteLine("UB:{0}, LB:{1}", UB, LB);
            #endregion
        }
        
    }
}
