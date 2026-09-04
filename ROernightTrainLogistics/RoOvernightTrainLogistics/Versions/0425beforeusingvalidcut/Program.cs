using System;
using System.Collections.Generic;
using System.Text;
using System.Linq;
using System.IO;
using System.Data;
using System.Diagnostics;
using System.Threading.Tasks;
using ILOG.Concert;
using ILOG.CPLEX;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.Distributions;
using System.Reflection;
//realize 1 DC distribution center

namespace RoOvernightTrainLogistics
{
    class Node
    {
        #region Node
        public string NodeName;        
        public int NodeID;  //external node number      
        public int NodeSeqNo;  //internal node number
        public int CorDcenterSeqNo;  //internal Dcenter number
        public double NodeDemand;//corresponding demand        
        public double DCstate;//1-as a candidate of DCs, otherwise 0
        public double AllocateLineNo;//line no for each node

        public List<Link> OutgoingLinkList = new List<Link>();//list of the arcs departure from the node
        public List<Link> IngoingLinkList = new List<Link>();//list of the arcs arrival the node

        #endregion
    }
    class DCenter
    {
        #region DCenter
        public string DCenterName;
        public int DCenterID;  //external node number      
        public int DCenterSeqNo;  //internal DCenter number
        public int CorNodeSeqNo;  //internal node number
        public double DCenterSupply;//corresponding Supply
        public double OpenCost;//corresponding cost if the node is open as a DC     
        public double NormalPrice;
        public double AllocateLineNo;//line no for each node

        public List<Link> OutgoingLinkList = new List<Link>();//list of the arcs departure from the node
        public List<Link> IngoingLinkList = new List<Link>();//list of the arcs arrival the node

        #endregion
    }
    class Link
    {
        #region Link

        public enum Type
        {
            DTD_LINK, DTC_LINK, CTC_LINK
        }
        public enum Direction
        {
            UP, DOWN
        }
        public int LinkID;// external arc id
        public double LinkTime;
        public double LinkCapacity;
        public double LinkCost;
        public int AllocateLineNo;
        public int CoreTSectionNo;
        public int CoreTSectionID;
        public Direction direction;

        public int LinkSeqNo;// internal arc id
        public int FromNodeID; // external from node id
        public int ToNodeID;// external to node id        
        public int InternalFromNodeSeqNo; // internal from node id,  
        public int InternalToNodeSeqNo;// internal to node id
        public string FromNodeName; // internal from node id,  
        public string ToNodeName;// internal to node id
        
        #endregion
    }
    class Section
    {
        #region Section

        public int SectionID;// external arc id
        public int SectionNo;// external arc id
        public double SectionCapacity;
        public List<Link> CorrespondLinks = new List<Link>();

        #endregion
    }
    class Path
    {
        #region Path
        public int PathNo;
        public double PathCost;        
        public double ReducedCost;
        public int xi_ij_coeff;
        public List<int> PathLinkSeqNoList = new List<int>();
        public List<int> PathNodeSeqNoList = new List<int>();
        public int StartNode;
        public int EndNode;

        #endregion
    }
    class Data
    {
        #region data
        public int DCSize;
        public int nodeSize;
        public int pathSize;
        public int linkSize;
        //public int sectionSize;
        //public List<Section> sectionlist;
        public List<Link> linklist;

        //After five consecutive iterations in which the LP bound does not improve, parameter λ is reset to 1
        public int Para_STB_IterLimit = 5;
        //After five more consecutive iterations with no LP bound improvement, parameter δ is reset to 0        
        //λ ∈ (0, 1] and δ ≥ 0 are initially set to 0.2 and 2ε(1e-5)
        public double Para_STB_Lambada = 0.2;
        public double Para_STB_Epsilon = 2 * 1e-5;
        public double Para_STB_Alpha = 0.9;
        public double Para_Scenario_fixthreshold = 0.2;
        public double Para_Location_fixthreshold = 0.2;
        public double para_scenario_prop = 0.02;//add 10% of columns obtained by CG to the second stage
        public bool Para_STB_strategy = true;//using strategy 1:  alpha * stab_y_sol[j] + (1 - alpha) * ast_y_sol[j]
        //otherwise, using strategy 2: λxi∗ +(1− λ) ˜xi + δ(1, . . . , 1), 

        public double[] f_j;
        public double[] c_j;
        public double[] B_i;

        //start and end node of a path
        public int[] s_r;
        public int[] e_r;
        
        public double[] d_r;
        public double[] h_i;
        public double MinCap;
        public double epsilon = 1e-6;//
        public double tolerance = 1e-5;//used for branch and price
        
        public int max_number_of_cols = 100;
        
        public int big_M;
        public double[] big_M_r;
      
        public int max_dstroyed_DCs;        
        public List<int>[] pathlist;
        //Technique 1: using techniques of benders decomposition
        //multiple cuts can be generated through two ways,
        //one for solving primal subproblem and enumerating all scenarios
        //the other for solving dual sub and using local branching
        
        public int maxMultiCuts = 3;
        public int localKappa = 1;//the parameter in local branching
        //Technique 2: using pareto strategy after identifying the worst scenario
        //at least a benders cut can be generated
        
        
        public bool multiCutStrategy = true;//true--using, false-unusing --STG1--
        public bool cuttingPlaneStrategy = true;//using cuttling plane at each node (except root node) --STG4--

        public bool paretoCutStrategy = false;//true--using, false-unusing --STG2-- for branch and benders cuts
        public bool stablilization_at_rootnode = false;//using stabilized strategy at root node--STG3--
        
        public bool validcut1 = false;
        public bool validcut2 = false;
        public bool subsetvalidcut = false;

        public bool solveCCGmasterIter = false;
                                                     
        public bool localBranchStrategy = false;//true--using, false-unusing --STG3-- for single level subproblem

        public bool BACforworstscenario = true;//true-using branch-and-cut for best scenario, solving MILP otherwise

        public bool solvingPrimalSub = true;////true--solving, false-unsolving
                                            ////
        public bool MultiplescenariosforCCG = false;

        public bool nogoodcuts = false;//whether or not to obtain dual values by solving primal problem
        
        public bool cplex_para_set = false;//set aggressive parameter to accelerate solution process   
        public double rootnodesolvetime = 10;

        public int[] pareto_xi;
        public double maxtransportcost;
        public double maxpurchasecost;
        public double maxopencost;

        public int[] pareto_y_sol;
        public int[] pareto_u_sol;

        public int[,] pathtoStartDC;
        public int[,] pathtoEndnode;
        public double LShapedbound;
        public double TL;

        public Data()
        {
            DCSize = Program.g_number_of_DCents;
            nodeSize = Program.g_number_of_nodes;
            pathSize = Program.g_number_of_paths;
            //sectionSize = Program.g_number_of_sections;
            max_dstroyed_DCs = 4;
            TL = 7200;

            linkSize = Program.g_number_of_links;

            pathtoStartDC = new int[pathSize, DCSize];
            pathtoEndnode = new int[pathSize, nodeSize];
            s_r = new int[pathSize];
            e_r = new int[pathSize];
            
            pathlist = new List<int>[pathSize];
            //sectionlist = Program.g_section_list;
            linklist = Program.g_link_list;

            for (int r = 0; r < pathSize; r++)
            {                
                s_r[r] = Program.g_node_list[Program.g_path_list[r].StartNode].CorDcenterSeqNo;
                e_r[r] = Program.g_path_list[r].EndNode;
                pathlist[r] = Program.g_path_list[r].PathLinkSeqNoList.ToList();
                pathtoStartDC[r, s_r[r]] = 1;
                pathtoEndnode[r, e_r[r]] = 1;
            }
            MinCap = int.MaxValue;
            for (int l = 0; l < linkSize; l++)
            {
                if(linklist[l].LinkCapacity < MinCap)
                {
                    MinCap = linklist[l].LinkCapacity;
                }
            }

            pareto_xi = new int[DCSize];

            c_j = new double[DCSize];
            
            f_j = new double[DCSize];

            d_r = new double[Program.g_number_of_paths];
            h_i = new double[nodeSize];

            for (int i = 0; i < nodeSize; i++)
            {
                h_i[i] = Program.g_node_list[i].NodeDemand;
            }
            maxtransportcost = 0;
            maxpurchasecost = 0;
            maxopencost = 0;

            for (int j = 0; j < DCSize; j++)
            {
                f_j[j] = Program.g_DCent_list[j].OpenCost;
                c_j[j] = Program.g_DCent_list[j].NormalPrice;

                if(maxopencost < f_j[j])
                {
                    maxopencost = f_j[j];
                }
                if(maxpurchasecost < c_j[j])
                {
                    maxpurchasecost = c_j[j];
                }
            }

            for (int l = 0; l < Program.g_number_of_paths; l++)
            {
                Path path = Program.g_path_list[l];

                d_r[l] = path.PathCost;

                if(maxtransportcost < d_r[l])
                {
                    maxtransportcost = d_r[l];
                }                
            }
            B_i = new double[nodeSize];
            for (int i = 0; i < nodeSize; i++)
            {
                B_i[i] = (maxtransportcost + maxpurchasecost/* + maxopencost*/) * h_i[i] * 2;
            }
            for (int i = 0; i < nodeSize; i++)
            {
                big_M += (int)Program.g_node_list[i].NodeDemand;
            }
            //they both are positive            
            big_M_r = new double[pathSize];
            for (int r = 0; r < pathSize; r++)
            {
                big_M_r[r] = (maxtransportcost + maxpurchasecost) * h_i[e_r[r]];
            }            
        }
        
        #endregion
    }
    class Solution
    {
        Data data = new Data();
        public double UB;
        public double LB;
        public double Gap;
        public double eps;

        int m_number_of_nodes = Program.g_number_of_nodes;
        int m_number_of_DCents = Program.g_number_of_DCents;
        int m_number_of_paths = Program.g_number_of_paths;

        public void write_solution(CCGMaster master)
        {
            #region write solution
            int[] y_solution = new int[m_number_of_DCents];
            double[] z_solution = new double[m_number_of_DCents];
            double[] t_solution = new double[m_number_of_DCents];
            double[,] x_ij_solution = new double[m_number_of_nodes, m_number_of_DCents];

            Program.solution_output.Write("DemandNodename" + ",");

            for (int i = 0; i < m_number_of_nodes; i++)
            {
                Program.solution_output.Write("{0}" + ",", Program.g_node_list[i].NodeName);

            }

            Program.solution_output.WriteLine();

            Program.solution_output.Write("DemandNodeID" + ",");

            for (int i = 0; i < m_number_of_nodes; i++)
            {
                Program.solution_output.Write("{0}" + ",", Program.g_node_list[i].NodeID);

            }
            Program.solution_output.WriteLine();
            Program.solution_output.Write("DCNodename" + ",");

            for (int j = 0; j < m_number_of_DCents; j++)
            {
                Program.solution_output.Write("{0}" + ",", Program.g_node_list[Program.g_DCent_list[j].CorNodeSeqNo].NodeName);

            }
            Program.solution_output.WriteLine();
            Program.solution_output.Write("DCNodeID" + ",");

            for (int j = 0; j < m_number_of_DCents; j++)
            {
                Program.solution_output.Write("{0}" + ",", Program.g_node_list[Program.g_DCent_list[j].CorNodeSeqNo].NodeID);

            }
            Program.solution_output.WriteLine();
            Program.solution_output.WriteLine();
            Program.solution_output.Write("Variables" + ",");

            for (int j = 0; j < m_number_of_DCents; j++)
            {
                Program.solution_output.Write("DC-{0}" + ",", Program.g_node_list[Program.g_DCent_list[j].CorNodeSeqNo].NodeName);

            }
            Program.solution_output.WriteLine();

            Program.solution_output.Write("y_sol" + ",");
            for (int j = 0; j < m_number_of_DCents; j++)
            {
                y_solution[j] = (int)master.model.GetValue(master.y_j[j]);
                Program.solution_output.Write("{0}" + ",", y_solution[j]);

            }
            Program.solution_output.WriteLine();

            //Program.solution_output.Write("z_sol" + ",");
            //for (int l = 0; l < master.z_l_j.Count; l++)
            //{
            //    for (int j = 0; j < m_number_of_DCents; j++)
            //    {
            //        z_solution[j] = master.model.GetValue(master.z_l_j[l][j]);
            //        Program.solution_output.Write("{0}" + ",", z_solution[j]);
            //    }
            //    Program.solution_output.WriteLine();
            //}
            //Program.solution_output.Write("s_sol" + ",");
            //int[] s_solution = new int[data.nodeSize];
            //for (int l = 0; l < master.s_l_i.Count; l++)
            //{
            //    for (int i = 0; i < m_number_of_nodes; i++)
            //    {
            //        s_solution[i] = (int) master.model.GetValue(master.s_l_i[l][i]);
            //        Program.solution_output.Write("{0}" + ",", s_solution[i]);
            //    }
            //    Program.solution_output.WriteLine();
            //}

            //Program.solution_output.WriteLine();

            //for (int l = 0; l < master.x_l_r.Count; l++)
            //{
            //    for (int r = 0; r < master.x_l_r[l].Length; r++)
            //    {
            //        x_ij_solution[data.e_r[r], data.s_r[r]] = master.model.GetValue(master.x_l_r[l][r]);
            //    }                

            //    for (int i = 0; i < m_number_of_nodes; i++)
            //    {                    
            //        Program.solution_output.Write("D-{0}" + ",", Program.g_node_list[i].NodeName);

            //        for (int j = 0; j < m_number_of_DCents; j++)
            //        {                        
            //            double nodedemand = 0;

            //            if (x_ij_solution[i, j] > 1e-6)
            //            {
            //                nodedemand = x_ij_solution[i, j] * Program.g_node_list[i].NodeDemand;
            //            }

            //            Program.solution_output.Write("{0}" + ",", nodedemand);
            //        }
            //        Program.solution_output.Write("{0}" + ",", Program.g_node_list[i].NodeDemand);
            //        Program.solution_output.WriteLine();
            //    }
            //    Program.solution_output.WriteLine();

            //    if (l == master.x_l_r.Count - 1)
            //    {                    

            //        double subobj = 0;
            //        for (int j = 0; j < data.DCSize; j++)
            //        {
            //            subobj += data.c_j[j] * z_solution[j];

            //        }
            //        for (int r = 0; r < data.pathSize; r++)
            //        {
            //            subobj += data.d_r[r] * data.h_i[data.e_r[r]] * x_ij_solution[data.e_r[r],data.s_r[r]];
            //        }
            //        for (int i = 0; i < data.nodeSize; i++)
            //        {
            //            subobj += data.B_i[i] * s_solution[i];
            //        }
            //        Console.WriteLine($"{master.model.GetValue(master.omega)},{subobj}");
            //    }
            //}
            #endregion
        }
        public void write_subproblemsolution(CCGSub subproblem)
        {
            #region write solution

            double[,] x_ij_solution = new double[data.nodeSize, data.DCSize];
            Cplex.Status solvstatus = subproblem.model.GetStatus();

            if(solvstatus == Cplex.Status.Optimal || solvstatus == Cplex.Status.Feasible)
            {
                Program.solution_output.Write("z_sol" + ",");
                for (int j = 0; j < m_number_of_DCents; j++)
                {
                    double z_val = subproblem.model.GetValue(subproblem.z_j[j]);
                    Program.solution_output.Write("{0}" + ",", z_val);
                }
                Program.solution_output.WriteLine();

                Program.solution_output.Write("x_sol" + ",");

                for (int j = 0; j < m_number_of_DCents; j++)
                {
                    Program.solution_output.Write("DC-{0}" + ",", Program.g_node_list[Program.g_DCent_list[j].CorNodeSeqNo].NodeID);

                }
                Program.solution_output.Write("NodeDemand");

                Program.solution_output.WriteLine();

                for (int r = 0; r < data.pathSize; r++)
                {
                    x_ij_solution[data.e_r[r], data.s_r[r]] = subproblem.model.GetValue(subproblem.x_r[r]);
                }

                for (int i = 0; i < m_number_of_nodes; i++)
                {
                    Program.solution_output.Write("Demand-{0}" + ",", Program.g_node_list[i].NodeName);

                    for (int j = 0; j < m_number_of_DCents; j++)
                    {
                        double nodedemand = 0;

                        if (x_ij_solution[i, j] > 1e-6)
                        {
                            nodedemand = x_ij_solution[i, j];
                        }

                        Program.solution_output.Write("{0}" + ",", nodedemand);
                    }


                    Program.solution_output.Write("{0}" + ",", Program.g_node_list[i].NodeDemand);
                    Program.solution_output.WriteLine();
                }
            }            
            #endregion
        }
        public void output_solution(BendersMaster master)
        {
            #region write solution
            int[] y_solution = new int[m_number_of_DCents];
            double[] z_solution = new double[m_number_of_DCents];            
            double[,] x_ij_solution = new double[m_number_of_nodes, m_number_of_DCents];

            Program.solution_output.Write("DemandNodename" + ",");

            for (int i = 0; i < m_number_of_nodes; i++)
            {
                Program.solution_output.Write("{0}" + ",", Program.g_node_list[i].NodeName);

            }

            Program.solution_output.WriteLine();

            Program.solution_output.Write("DemandNodeID" + ",");

            for (int i = 0; i < m_number_of_nodes; i++)
            {
                Program.solution_output.Write("{0}" + ",", Program.g_node_list[i].NodeID);

            }
            Program.solution_output.WriteLine();
            Program.solution_output.Write("DCNodename" + ",");

            for (int j = 0; j < m_number_of_DCents; j++)
            {
                Program.solution_output.Write("{0}" + ",", Program.g_node_list[Program.g_DCent_list[j].CorNodeSeqNo].NodeName);

            }
            Program.solution_output.WriteLine();
            Program.solution_output.Write("DCNodeID" + ",");

            for (int j = 0; j < m_number_of_DCents; j++)
            {
                Program.solution_output.Write("{0}" + ",", Program.g_node_list[Program.g_DCent_list[j].CorNodeSeqNo].NodeID);

            }
            Program.solution_output.WriteLine();
            Program.solution_output.WriteLine();
            Program.solution_output.Write("Variables" + ",");

            for (int j = 0; j < m_number_of_DCents; j++)
            {
                Program.solution_output.Write("DC-{0}" + ",", Program.g_node_list[Program.g_DCent_list[j].CorNodeSeqNo].NodeName);

            }
            Program.solution_output.WriteLine();

            Program.solution_output.Write("y_sol" + ",");

            for (int j = 0; j < m_number_of_DCents; j++)
            {
                var yval = master.model.GetValue(master.y_j[j]);

                if (yval > 0.5)
                {
                    y_solution[j] = 1;
                }

                Program.solution_output.Write("{0}" + ",", y_solution[j]);

            }
            Program.solution_output.WriteLine();

            //Program.solution_output.Write("r_sol" + ",");
            //for (int j = 0; j < m_number_of_DCents; j++)
            //{
            //    double r_solution = subproblem.model.GetValue(subproblem.r_j[j]);
            //    Program.solution_output.Write("{0}" + ",", r_solution);
            //}
            //Program.solution_output.WriteLine();


            //Program.solution_output.Write("t_sol" + ",");
            //for (int j = 0; j < m_number_of_DCents; j++)
            //{
            //    t_solution[j] = subproblem.model.GetValue(subproblem.t_j[j]);

            //    double nodedemand = 0;

            //    if (t_solution[j] > 1e-6)
            //    {
            //        nodedemand = t_solution[j];
            //    }

            //    Program.solution_output.Write("{0}" + ",", nodedemand);

            //}

            //Program.solution_output.WriteLine();

            //Program.solution_output.Write("x_sol" + ",");

            //for (int j = 0; j < m_number_of_DCents; j++)
            //{
            //    Program.solution_output.Write("DC-{0}" + ",", Program.g_node_list[Program.g_DCent_list[j].CorNodeSeqNo].NodeID);

            //}
            //Program.solution_output.Write("NodeDemand");

            //Program.solution_output.WriteLine();

            //for (int i = 0; i < m_number_of_nodes; i++)
            //{
            //    Program.solution_output.Write("Demand-{0}" + ",", Program.g_node_list[i].NodeName);

            //    for (int j = 0; j < m_number_of_DCents; j++)
            //    {
            //        x_ij_solution[i, j] = subproblem.model.GetValue(subproblem.x_ij[i, j]);
            //        double nodedemand = 0;

            //        if (x_ij_solution[i, j] > 1e-6)
            //        {
            //            nodedemand = x_ij_solution[i, j] * Program.g_node_list[i].NodeDemand;
            //        }

            //        Program.solution_output.Write("{0}" + ",", nodedemand);
            //    }


            //    Program.solution_output.Write("{0}" + ",", Program.g_node_list[i].NodeDemand);
            //    Program.solution_output.WriteLine();
            //}
            #endregion
        }
        public void outputSubproSolution(BendersMaster master)
        {
            #region write solution
            int[] y_solution = new int[m_number_of_DCents];
            int[] xi_sol = new int[m_number_of_DCents];
            Data data = new Data();

            for (int j = 0; j < m_number_of_DCents; j++)
            {
                y_solution[j] = (int)master.model.GetValue(master.y_j[j]);                
            }

            CACGBD cobdgen = new CACGBD();

            SCNRMaster submaster = cobdgen.BranchandCutForWorstCsenario(y_solution);

            Cplex.Status cur_status = submaster.model.GetStatus();

            if (!submaster.feasiblestatus)
            {
                xi_sol = submaster.worstScenarioSolution.ToArray();
            }
            else
            {
                
                for (int l = 0; l < data.DCSize; l++)
                {
                    var u_val = submaster.model.GetValue(submaster.u_j[l]);

                    if (u_val > 0.5)
                    {
                        xi_sol[l] = 1;
                    }
                }
            }

            CCGSub subproblem = new CCGSub();

            subproblem.GenCCGSubproblem(y_solution, xi_sol);
            subproblem.model.Solve();
            
            double[,] x_ij_solution = new double[m_number_of_nodes, m_number_of_DCents];

            Program.solution_output.Write("z_sol" + ",");
            for (int j = 0; j < m_number_of_DCents; j++)
            {
                double z_solution = subproblem.model.GetValue(subproblem.z_j[j]);
                Program.solution_output.Write("{0}" + ",", z_solution);
            }
            Program.solution_output.WriteLine();

            Program.solution_output.Write("x_sol" + ",");

            for (int j = 0; j < m_number_of_DCents; j++)
            {
                Program.solution_output.Write("DC-{0}" + ",", Program.g_node_list[Program.g_DCent_list[j].CorNodeSeqNo].NodeID);

            }
            Program.solution_output.Write("NodeDemand");

            Program.solution_output.WriteLine();

            for (int i = 0; i < m_number_of_nodes; i++)
            {
                Program.solution_output.Write("Demand-{0}" + ",", Program.g_node_list[i].NodeName);

                for (int j = 0; j < m_number_of_DCents; j++)
                {
                    for (int r = 0; r < m_number_of_paths; r++)
                    {
                        if (data.e_r[r] == i && data.s_r[r] == j)
                        {
                            x_ij_solution[i, j] = subproblem.model.GetValue(subproblem.x_r[r]);
                            double nodedemand = 0;

                            if (x_ij_solution[i, j] > 1e-6)
                            {
                                nodedemand = x_ij_solution[i, j] /** Program.g_node_list[i].NodeDemand*/;
                            }

                            Program.solution_output.Write("{0}" + ",", nodedemand);
                        }
                    }
                    
                }


                Program.solution_output.Write("{0}" + ",", Program.g_node_list[i].NodeDemand);
                Program.solution_output.WriteLine();
            }
            #endregion
        }
    }    
    
    class Program
    {
        #region global parameters
        public static int g_number_of_nodes;
        public static int g_number_of_links;
        public static int g_number_of_DCents;
        public static int g_number_of_paths;
        public static int g_number_of_sections;

        public static int CCG_debug = 0;
        public static int transportmode = 0;//0-road transportation,1-helicopter

        public static List<Node> g_node_list;
        public static List<Link> g_link_list;
        public static List<Section> g_section_list;
        public static List<DCenter> g_DCent_list;
        public static List<Path> g_path_list;
        
        public static Dictionary<int, int> g_node_id_to_internal_node_seq_no_dic;
        public static Dictionary<int, int> g_internal_node_seq_no_to_node_id_dic;

        public static Dictionary<int, int> g_link_id_to_internal_link_seq_no_dic;
        public static Dictionary<int, int> g_internal_link_seq_no_to_link_id_dic;
        public static Dictionary<int, int> g_section_id_to_internal_section_seq_no_dic;

        //.csv documents used for saving information of nodes, links, and the solution time so on
        static string s1 = @"..\..\..\..\TestLog\g_pFileOutputLog.csv";
        public static StreamWriter g_pFileOutputLog = new StreamWriter(s1, false, Encoding.ASCII);
        static string s5 = @"..\..\..\..\TestLog\g_parameteroutput.csv";
        public static StreamWriter g_parameteroutput = new StreamWriter(s5, false, Encoding.ASCII);
        static string s6 = @"..\..\..\..\TestLog\solution_output.csv";
        public static StreamWriter solution_output = new StreamWriter(s6, false, Encoding.ASCII);
        static string s9 = @"..\..\..\..\TestLog\SolutionIteration.csv";
        public static StreamWriter SolutionIteration = new StreamWriter(s9, false, Encoding.ASCII);
        static string s10 = @"..\..\..\..\TestLog\g_comparativeIndicators.csv";
        public static StreamWriter g_comparativeIndicators = new StreamWriter(s10, false, Encoding.ASCII);
        static string s11 = @"..\..\..\..\TestLog\g_CCGiteration.csv";
        public static StreamWriter g_CCGiteration = new StreamWriter(s11, false, Encoding.ASCII);
        static string s12 = @"..\..\..\..\TestLog\g_CPLEXResults.csv";
        public static StreamWriter g_CPLEXResults = new StreamWriter(s12, false, Encoding.ASCII);

        static string s2 = @"..\..\..\..\TestLog\RunOutput.txt";
        public static TextWriter TWoutput = File.CreateText(s2);
        static string s3 = @"..\..\..\..\TestLog\CCG_debugfile.txt";
        public static TextWriter CCG_debugfile = File.CreateText(s3);
        static string s4 = @"..\..\..\..\TestLog\CCG_resultsfile.txt";
        public static TextWriter CCG_resultsfile = File.CreateText(s4);
        static string s7 = @"..\..\..\..\TestLog\BDD_debugfile.txt";
        public static TextWriter BDD_debugfile = File.CreateText(s7);
        static string s8 = @"..\..\..\..\TestLog\Solvinglog.txt";
        public static TextWriter Solvinglog = File.CreateText(s8);

        #endregion

        /// <summary>
        /// Step one, reading data
        /// </summary>
        public static void g_ReadInputData()
        {
            #region reading

            #region variables, list and dics

            //initialize global variables            
            g_number_of_nodes = 0;
            g_number_of_links = 0;
            g_number_of_DCents = 0;
            g_number_of_paths = 0;
            g_number_of_sections = 0;

            //initialize dic
            g_node_id_to_internal_node_seq_no_dic = new Dictionary<int, int>();  // dictionary, map node number to internal node sequence no. 
            g_internal_node_seq_no_to_node_id_dic = new Dictionary<int, int>();  // dictionary, map internal node sequence no to node number.

            g_link_id_to_internal_link_seq_no_dic = new Dictionary<int, int>();
            g_internal_link_seq_no_to_link_id_dic = new Dictionary<int, int>();

            g_section_id_to_internal_section_seq_no_dic = new Dictionary<int, int>();

            //initialize list
            g_node_list = new List<Node>();//saving all parameters from nodes
            g_link_list = new List<Link>();//saving all parameters from links
            g_DCent_list = new List<DCenter>();
            g_path_list = new List<Path>();
            g_section_list = new List<Section>();

            #endregion

            string[] fileName = { "Input_Nodes", "Input_Links", "Input_Paths"};

            for (int s = 0; s < fileName.Length; s++)
            {
                string filePath = @"..\..\..\..\Dataset\";
                filePath += fileName[s];
                filePath += ".csv";

                FileStream fs = new FileStream(filePath, FileMode.Open, FileAccess.Read);
                StreamReader sr = new StreamReader(fs, System.Text.Encoding.Default);

                if (s == 0)//代表输入点的信息
                {

                    //记录每次读取的一行记录
                    string strLine = "";

                    //记录每行记录中的各字段内容
                    string[] aryLine;

                    //标示列数
                    int columnCount = 0;

                    //标示是否是读取的第一行
                    bool IsFirst = true;

                    //逐行读取CSV中的数据
                    while ((strLine = sr.ReadLine()) != null)
                    {
                        string name;
                        int node_id, DCstate, allocateno;
                        double demand, cost, nprice, reprice, emprice;

                        aryLine = strLine.Split(',');

                        if (IsFirst == true)
                        {
                            IsFirst = false;
                            columnCount = aryLine.Length;
                            continue;

                        }
                        else
                        {
                            node_id = int.Parse(aryLine[0]); name = aryLine[1]; demand = double.Parse(aryLine[2]); DCstate = int.Parse(aryLine[3]);
                            cost = double.Parse(aryLine[4]); nprice = double.Parse(aryLine[5]); allocateno = int.Parse(aryLine[6]);

                            g_internal_node_seq_no_to_node_id_dic.Add(g_number_of_nodes, node_id);
                            g_node_id_to_internal_node_seq_no_dic.Add(node_id, g_number_of_nodes);


                            Node node = new Node();//实例化一个node
                            node.NodeName = name;
                            node.NodeID = node_id;//虚拟小区编号
                            node.NodeDemand = demand;//小区需求                            
                            node.NodeSeqNo = g_number_of_nodes;
                            node.DCstate = DCstate;
                            node.AllocateLineNo = allocateno;

                            if (DCstate == 1)
                            {
                                DCenter dCenter = new DCenter();

                                dCenter.DCenterName = name;
                                dCenter.DCenterID = node_id;//虚拟小区编号                                                       
                                dCenter.DCenterSeqNo = g_number_of_DCents;
                                dCenter.CorNodeSeqNo = g_number_of_nodes;
                                dCenter.OpenCost = cost;
                                dCenter.NormalPrice = nprice;
                                dCenter.AllocateLineNo = allocateno;

                                node.CorDcenterSeqNo = g_number_of_DCents;

                                g_DCent_list.Add(dCenter);
                                g_number_of_DCents++;
                            }
                            else
                            {
                                node.CorDcenterSeqNo = -1;
                            }

                            g_number_of_nodes++;
                            g_node_list.Add(node);

                            if (g_number_of_nodes % 10 == 0)
                                Console.WriteLine("reading" + " " + g_number_of_nodes + " " + "nodes");
                        }
                    }

                    Console.WriteLine("The number of nodes =" + " " + g_number_of_nodes + " " + "nodes");
                    g_pFileOutputLog.WriteLine("The number of nodes" + "," + g_number_of_nodes);

                    Console.WriteLine("The number of dCenters =" + " " + g_number_of_DCents + " " + "dCenters");
                    g_pFileOutputLog.WriteLine("The number of dCenters" + "," + g_number_of_DCents);
                }
                else if (s == 1)
                {
                    //记录每次读取的一行记录
                    string strLine = "";

                    //记录每行记录中的各字段内容
                    string[] aryLine;

                    //标示列数
                    int columnCount = 0;

                    //标示是否是读取的第一行
                    bool IsFirst = true;

                    //逐行读取CSV中的数据
                    while ((strLine = sr.ReadLine()) != null)
                    {
                        string from_node_name, to_node_name, direct;
                        int link_id, link_time, from_node_id, to_node_id, allocatelinno, sectionID;
                        double link_cost,capacity;

                        aryLine = strLine.Split(',');

                        if (IsFirst == true)
                        {
                            IsFirst = false;
                            columnCount = aryLine.Length;
                            continue;

                        }
                        else
                        {
                            link_id = int.Parse(aryLine[0]); from_node_id = int.Parse(aryLine[1]); to_node_id = int.Parse(aryLine[2]); from_node_name = aryLine[3];
                            to_node_name = aryLine[4]; link_time = int.Parse(aryLine[5]); link_cost = double.Parse(aryLine[6]); allocatelinno = int.Parse(aryLine[7]);
                            sectionID = int.Parse(aryLine[8]); direct = aryLine[9]; capacity = double.Parse(aryLine[10]);
                            Link link = new Link();//实例化一个Link

                            link.LinkID = link_id;//路段的id对应弧的起点
                            link.FromNodeID = from_node_id;//路段的起点node_id
                            link.ToNodeID = to_node_id;//路段的终点node_id                            

                            int internal_from_node_seq_no = g_node_id_to_internal_node_seq_no_dic[link.FromNodeID];  // 路段的起点node_id对应的内部编号
                            int internal_to_node_seq_no = g_node_id_to_internal_node_seq_no_dic[link.ToNodeID];// 路段的终点node_id对应的内部编号 

                            link.InternalFromNodeSeqNo = internal_from_node_seq_no;//路段的内部起点
                            link.InternalToNodeSeqNo = internal_to_node_seq_no;//路段的内部终点                            

                            link.FromNodeName = from_node_name;//路段起点名称
                            link.ToNodeName = to_node_name;//路段终点名称

                            link.LinkSeqNo = g_number_of_links;//路段编号

                            link.LinkTime = link_time;//路段时间
                            link.LinkCost = link_cost; ;//运输费用
                            link.AllocateLineNo = allocatelinno;//所属线路
                            link.CoreTSectionID = sectionID;
                            link.LinkCapacity = capacity;

                            //该点的出发弧集合
                            g_node_list[internal_from_node_seq_no].OutgoingLinkList.Add(link);//把路段加到从该点出发的路段列表
                            g_node_list[internal_to_node_seq_no].IngoingLinkList.Add(link);//把路段加入到达该点的路段列表

                            if (g_node_list[internal_from_node_seq_no].CorDcenterSeqNo != -1)
                            {
                                g_DCent_list[g_node_list[internal_from_node_seq_no].CorDcenterSeqNo].OutgoingLinkList.Add(link);//把路段加到从该点出发的路段列表
                            }
                            if (g_node_list[internal_to_node_seq_no].CorDcenterSeqNo != -1)
                            {
                                g_DCent_list[g_node_list[internal_to_node_seq_no].CorDcenterSeqNo].IngoingLinkList.Add(link);//把路段加入到达该点的路段列表
                            }

                            if(direct == "UP")
                            {
                                link.direction = Link.Direction.UP;
                                Section section = new Section();

                                section.SectionID = sectionID;
                                section.SectionNo = g_number_of_sections;
                                section.CorrespondLinks.Add(link);
                                section.SectionCapacity = capacity;

                                link.CoreTSectionNo = g_number_of_sections;//所属线路
                                g_section_id_to_internal_section_seq_no_dic.Add(sectionID, g_number_of_sections);

                                g_section_list.Add(section);
                                g_number_of_sections++;
                            }
                            else
                            {
                                link.direction = Link.Direction.DOWN;
                            }

                            g_link_id_to_internal_link_seq_no_dic.Add(link.LinkID, link.LinkSeqNo);//外部和内部关联
                            g_internal_link_seq_no_to_link_id_dic.Add(link.LinkSeqNo, link.LinkID);//内部和外部关联 


                            g_number_of_links++;
                            g_link_list.Add(link);

                            if (g_number_of_links % 50 == 0)
                                Console.WriteLine("reading" + " " + g_number_of_links + " " + "links");
                        }
                    }

                    Console.WriteLine("The number of links =" + " " + g_number_of_links + " " + "links");
                    g_pFileOutputLog.WriteLine("The number of links" + "," + g_number_of_links);

                    for (int l = 0; l < g_number_of_links; l++)
                    {
                        Link link = g_link_list[l];
                        int sectionID = link.CoreTSectionID;

                        if(link.direction == Link.Direction.DOWN)
                        {
                            int sectionno = g_section_id_to_internal_section_seq_no_dic[sectionID];
                            Section section = g_section_list[sectionno];
                            section.CorrespondLinks.Add(link);
                            link.CoreTSectionNo = sectionno;
                        }
                    }
                }
                else if (s == 2)
                {

                    List<List<double>> temp_mtr = new List<List<double>>();
                    //记录每次读取的一行记录
                    string strLine = "";

                    //记录每行记录中的各字段内容
                    string[] aryLine;

                    //标示列数
                    int columnCount = 0;

                    //标示是否是读取的第一行
                    bool IsFirst = true;

                    //逐行读取CSV中的数据
                    while ((strLine = sr.ReadLine()) != null)
                    {
                        aryLine = strLine.Split(',');

                        int path_id, from_dcentid, to_communityid;
                        double cost;

                        if (IsFirst == true)
                        {
                            IsFirst = false;
                            columnCount = aryLine.Length;
                            continue;

                        }
                        else
                        {
                            path_id = int.Parse(aryLine[0]); from_dcentid = int.Parse(aryLine[1]); to_communityid = int.Parse(aryLine[2]); cost = double.Parse(aryLine[3]);

                            Path path = new Path();
                            path.PathNo = path_id;
                            path.StartNode = g_node_id_to_internal_node_seq_no_dic[from_dcentid];
                            path.EndNode = g_node_id_to_internal_node_seq_no_dic[to_communityid];
                            path.PathCost = cost;

                            if (from_dcentid != to_communityid)
                            {
                                for (int p = 4; p < aryLine.Length; p++)
                                {
                                    if(aryLine[p] != "")
                                    {
                                        path.PathLinkSeqNoList.Add(g_link_id_to_internal_link_seq_no_dic[int.Parse(aryLine[p])]);
                                    }                                    
                                }
                            }

                            g_number_of_paths++;
                            g_path_list.Add(path);

                        }
                        if (g_number_of_paths % 50 == 0)
                            Console.WriteLine("reading" + " " + g_number_of_paths + " " + "paths");
                    }
                    Console.WriteLine("The number of paths =" + " " + g_number_of_paths + " " + "paths");
                    g_pFileOutputLog.WriteLine("The number of paths" + "," + g_number_of_paths);
                }                

            }
            #endregion
        }        
        /// <summary>
        /// 第三步Column generation and ADMM下界求解
        /// </summary>
        /// <returns></returns>
        public static void g_RoMaxFacilityLocation_optimization()
        {
            #region 调用 CCG & CG 求解
            Stopwatch optimize_MILP = new Stopwatch();
            Stopwatch optimize_CCGwBAC = new Stopwatch();
            Stopwatch optimize_BABC = new Stopwatch();
            Stopwatch optimize_ParaCCG = new Stopwatch();

            float unit = 1000F;
            float total_CCGBAC = 0;
            float total_MILP = 0;
            float total_BABC = 0;
            float total_ParaCCG = 0;

            Data dt = new Data();
            //stage 1 : column and constraint generation optimization
            int ParaCCG_debug_flag = 2;
            if (ParaCCG_debug_flag == 1)
            {
                Solution solution = new Solution();

                optimize_ParaCCG.Start();
                CandCG candCG = new CandCG();
                candCG.ParallelCCG();

                optimize_ParaCCG.Stop();
                total_ParaCCG = optimize_ParaCCG.ElapsedMilliseconds;
            }

            int MILPCCG_debug_flag = 2;
            if (MILPCCG_debug_flag == 1)
            {
                Solution solution = new Solution();
                
                optimize_MILP.Start();
                CandCG candCG = new CandCG();
                candCG.CCGMILP();               

                optimize_MILP.Stop();
                total_MILP = optimize_MILP.ElapsedMilliseconds;
            }

            int CCGBAC_debug_flag = 2;
            if (CCGBAC_debug_flag == 1)
            {
                Solution solution = new Solution(); 

                optimize_CCGwBAC.Start();
                CandCG candCG = new CandCG();
                candCG.CCGwBAC();

                optimize_CCGwBAC.Stop();
                total_CCGBAC = optimize_CCGwBAC.ElapsedMilliseconds;
            }

            //stage 2 : column and Benders cut optimization
            int CGBD_debug_flag = 1;
            if (CGBD_debug_flag == 1)
            {
                optimize_BABC.Start();

                Solution solution = new Solution();
                CACGBD colBDGen = new CACGBD();

                BendersMaster Bendersmaster = colBDGen.BranchandCut();
                optimize_BABC.Stop();
                total_BABC = optimize_BABC.ElapsedMilliseconds;

                solution.output_solution(Bendersmaster); solution.outputSubproSolution(Bendersmaster);
                Console.WriteLine($"Optimal master objective cost: {Bendersmaster.model.GetObjValue()}");
                Console.WriteLine($"Best lower bound: {Bendersmaster.model.GetBestObjValue()}");
            }

            Console.WriteLine("CPU Running Time of total_ParaCCG = {0} seconds", total_ParaCCG / unit);
            Console.WriteLine("CPU Running Time of total_MILPCCG = {0} seconds", total_MILP / unit);
            Console.WriteLine("CPU Running Time of total_CCGBAC = {0} seconds", total_CCGBAC / unit);
            Console.WriteLine("CPU Running Time of BABC = {0} seconds", total_BABC / unit);

            g_pFileOutputLog.Write("CPU Running Time of total_ParaCCG =,{0}, seconds\r\n", total_ParaCCG / unit);
            g_pFileOutputLog.Write("CPU Running Time of total_MILPCCG =,{0}, seconds\r\n", total_MILP / unit);
            g_pFileOutputLog.Write("CPU Running Time of total_CCGBAC =,{0}, seconds\r\n", total_CCGBAC / unit);
            g_pFileOutputLog.Write("CPU Running Time of BABC =,{0}, seconds\r\n", total_BABC / unit);

            g_pFileOutputLog.WriteLine("STG1-multicuts, STG2-Paretocut, STG3-Localbranch,STG4-Stabilization");
            g_pFileOutputLog.WriteLine($"{dt.multiCutStrategy}, {dt.paretoCutStrategy}, {dt.localBranchStrategy},{dt.stablilization_at_rootnode}");

            g_pFileOutputLog.WriteLine("cplexparaset, solvingprimalsub, nogoodcuts,BACforworstscenario,MultiplescenariosforCCG");
            g_pFileOutputLog.WriteLine($"{dt.cplex_para_set}, {dt.solvingPrimalSub}, {dt.nogoodcuts},{dt.BACforworstscenario},{dt.MultiplescenariosforCCG}");

            #endregion
        }
        static void Main(string[] args)
        {
            // step 1: read input data of network and demand agent
            g_ReadInputData();
            
            //generateMatrix();            

            g_RoMaxFacilityLocation_optimization();

            g_parameteroutput.WriteLine();

            g_pFileOutputLog.Flush();
            g_pFileOutputLog.Close();
            CCG_debugfile.Flush();
            CCG_debugfile.Close();
            CCG_resultsfile.Flush();
            CCG_resultsfile.Close();
            g_parameteroutput.Flush();
            g_parameteroutput.Close();
            solution_output.Flush();
            solution_output.Close();

            BDD_debugfile.Flush();
            BDD_debugfile.Close();
            Solvinglog.Flush();
            Solvinglog.Close();
            SolutionIteration.Flush();
            SolutionIteration.Close();

            g_comparativeIndicators.Flush();
            g_comparativeIndicators.Close();

            g_CPLEXResults.Flush();
            g_CPLEXResults.Close();

            g_CCGiteration.Flush();
            g_CCGiteration.Close();

            Console.WriteLine("End of Optimization ");            
            Console.WriteLine("done.");
            Console.ReadKey();
        }
    }
}
