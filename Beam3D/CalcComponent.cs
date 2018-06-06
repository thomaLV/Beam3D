using System;
using System.Collections.Generic;

using Grasshopper.Kernel;
using Rhino.Geometry;
using System.Drawing;
using Grasshopper.GUI.Canvas;
using System.Windows.Forms;
using Grasshopper.GUI;

using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;
using System.Diagnostics;

namespace Beam3D
{
    public class CalcComponent : GH_Component
    {
        public CalcComponent()
          : base("BeamCalculation", "BeamC",
              "Description",
              "Koala", "3D Beam")
        {
        }

        //Initialize moments
        static bool startCalc = false;
        //static bool startTest = false;

        //Method to allow c hanging of variables via GUI (see Component Visual)
        public static void setStart(string s, bool i)
        {
            if (s == "Run")
            {
                startCalc = i;
            }
            //if (s == "Run Test")
            //{
            //    startTest = i;
            //}
        }

        public override void CreateAttributes()
        {
            m_attributes = new Attributes_Custom(this);
        }

        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddLineParameter("Lines", "LNS", "Geometry, in form of Lines)", GH_ParamAccess.list);
            pManager.AddTextParameter("Boundary Conditions", "BDC", "Boundary Conditions in form x,y,z,vx,vy,vz,rx,ry,rz", GH_ParamAccess.list);
            pManager.AddTextParameter("Material properties", "Mat", "Material Properties: E, A, Iy, Iz, v, alpha (rotation about x)", GH_ParamAccess.item, "200000,3600,4920000,4920000, 0.3, 0");
            pManager.AddTextParameter("PointLoads", "PL", "Load given as Vector [N]", GH_ParamAccess.list);
            pManager.AddTextParameter("PointMoment", "PM", "Moment set in a point in [Nm]", GH_ParamAccess.list, "");
            pManager.AddIntegerParameter("Elements", "n", "Number of elements", GH_ParamAccess.item, 1);
        }

        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            //pManager.AddGenericParameter("Deformations", "Def", "Tree of Deformations", GH_ParamAccess.tree);
            pManager.AddNumberParameter("Deformations", "Def", "Tree of Deformations", GH_ParamAccess.list);
            pManager.AddNumberParameter("Reaction Forces", "R", "Reaction Forces", GH_ParamAccess.list);
            pManager.AddNumberParameter("Applied Loads", "A", "Applied Loads", GH_ParamAccess.list);
            //pManager.AddGenericParameter("Element stresses", "Strs", "The Stress in each element", GH_ParamAccess.tree);
            //pManager.AddGenericParameter("Element strains", "Strn", "The Strain in each element", GH_ParamAccess.tree);
            pManager.AddNumberParameter("Element stresses", "Strs", "The Stress in each element", GH_ParamAccess.list);
            pManager.AddNumberParameter("Element strains", "Strn", "The Strain in each element", GH_ParamAccess.list);
            pManager.AddGenericParameter("Matrix Deformations", "DM", "Deformation Matrix for def. component", GH_ParamAccess.item);
            pManager.AddPointParameter("New Base Points", "NBP", "Nodal points of sub elements", GH_ParamAccess.list);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
           #region Fetch input
            //Expected inputs
            List<Line> geometry = new List<Line>();         //Initial Geometry of lines
            List<string> bdctxt = new List<string>();       //Boundary conditions in string format
            List<string> loadtxt = new List<string>();      //loads in string format
            List<string> momenttxt = new List<string>();    //Moments in string format
            string mattxt = "";
            int n = 1;


            //Set expected inputs from Indata
            if (!DA.GetDataList(0, geometry)) return;       //sets geometry
            if (!DA.GetDataList(1, bdctxt)) return;         //sets boundary conditions as string
            if (!DA.GetData(2, ref mattxt)) return;         //sets material properties as string
            if (!DA.GetDataList(3, loadtxt)) return;        //sets load as string
            if (!DA.GetDataList(4, momenttxt)) return;      //sets moment as string
            if (!DA.GetData(5, ref n)) return;              //sets number of elements
            #endregion

            //Interpret and set material parameters
            double E;       //Material Young's modulus, initial value 200 000 [MPa]
            double A;       //Area for each element in same order as geometry, initial value CFS100x100 3600 [mm^2]
            double Iy;      //Moment of inertia about local y axis, initial value 4.92E6 [mm^4]
            double Iz;      //Moment of inertia about local z axis, initial value 4.92E6 [mm^4]
            double J;       //Polar moment of inertia
            double G;       //Shear modulus, initial value 79300 [mm^4]
            double v;       //Poisson's ratio, initial value 0.3
            double alpha;


            SetMaterial(mattxt, out E, out A, out Iy, out Iz, out J, out G, out v, out alpha);

            #region Prepares geometry, boundary conditions and loads for calculation
            //List all nodes (every node only once), numbering them according to list index
            List<Point3d> points = CreatePointList(geometry);


            //Interpret the BDC inputs (text) and create list of boundary condition (1/0 = free/clamped) for each dof.
            Vector<double> bdc_value = CreateBDCList(bdctxt, points);


            //Interpreting input load (text) and creating load list (do uble)
            Vector<double> load = CreateLoadList(loadtxt, momenttxt, points);
            #endregion

            Matrix<double> def_shape, glob_strain, glob_stress;
            Vector<double> reactions;
            List<Point3d> oldXYZ;
            
            List<Curve> defGeometry = new List<Curve>();    //output deformed geometry


            if (startCalc)
            {
                #region Create global and reduced stiffness matrix
                //Create global stiffness matrix
                Matrix<double> K_tot = GlobalStiffnessMatrix(geometry, points, E, A, Iy, Iz, J, G, alpha);


                //Create reduced K-matrix and reduced load list (removed free dofs)
                Matrix<double> KGr;
                Vector<double> load_red;
                ReducedGlobalStiffnessMatrix(bdc_value, K_tot, load, out KGr, out load_red);
                #endregion

                #region Solver Performance Test
                //if (startTest)
                //{
                //    string output_time = "";
                //    string performanceResult = "=================START OF TEST=================" + Environment.NewLine;
                //    performanceResult += "Number of lines: " + geometry.Count.ToString() + Environment.NewLine;
                //    string tester = "";

                //    //checking for error with writing to file (skip test if unable to write)
                //    try
                //    {
                //        using (System.IO.StreamWriter file =
                //        new System.IO.StreamWriter(@"solverTest.txt", true))
                //        {
                //            file.WriteLine(tester);
                //        }
                //    }
                //    //create new file if solverTest.txt does not exist
                //    catch (System.IO.DirectoryNotFoundException)
                //    {
                //        System.IO.File.WriteAllText(@"solverTest.txt", tester);
                //    }
                //    //other exception (no write access?)
                //    catch (Exception)
                //    {
                //        AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Write to file error! Create file solverTest.txt in Koala folder and/or Rhinoceros 5.exe folder");
                //        return;
                //    }

                //    int decimals = 2;
                //    CheckSolvers(KGr, load_red, decimals, out output_time);

                //    performanceResult += output_time;

                //    //append result to txt-file (at buildpath)
                //    try
                //    {
                //        using (System.IO.StreamWriter file =
                //        new System.IO.StreamWriter(@"solverTest.txt", true))
                //        {
                //            file.WriteLine(performanceResult);
                //        }
                //    }
                //    //create new file if solverTest.txt does not exist
                //    catch (System.IO.DirectoryNotFoundException)
                //    {
                //        System.IO.File.WriteAllText(@"\solverTest.txt", output_time);
                //    }

                //}
                #endregion

                #region Calculate deformations, reaction forces and internal strains and stresses
                //Calculate deformations
                Vector<double> def_red = KGr.Cholesky().Solve(load_red);


                //Add the clamped dofs (= 0) to the deformations list
                Vector<double> def_tot = RestoreTotalDeformationVector(def_red, bdc_value);


                //Calculate the reaction forces from the deformations
                reactions = K_tot.Multiply(def_tot);
                reactions -= load; //method for separating reactions and applied loads
                reactions.CoerceZero(1e-10);


                //Interpolate deformations using shape functions
                double y = 50;

                var z = y;
                InterpolateDeformations(def_tot, points, geometry, n, z, y, alpha, out def_shape, out oldXYZ, out glob_strain);


                //Calculate stresses
                glob_stress = E * glob_strain;
                #endregion

            }
            else
            {
                #region Set outputs to zero
                reactions = Vector<double>.Build.Dense(points.Count * 6);
                def_shape = Matrix<double>.Build.Dense(geometry.Count, 6 * (n + 1));
                glob_strain = def_shape;
                glob_stress = def_shape;

                oldXYZ = new List<Point3d>();
                #endregion
            }

            #region Format output
            //Grasshopper.DataTree<double> def_shape_nested = ConvertToNestedList(def_shape);
            //Grasshopper.DataTree<double> strain_nested = ConvertToNestedList(glob_strain);
            //Grasshopper.DataTree<double> stresses_nested = ConvertToNestedList(glob_stress);

            double[,] def_m = new double[def_shape.RowCount, def_shape.ColumnCount];
            for (int i = 0; i < def_shape.RowCount; i++)
            {
                for (int j = 0; j < def_shape.ColumnCount; j++)
                {
                    def_m[i, j] = def_shape[i, j];
                }
            }

            double[] def = new double[def_shape.RowCount * def_shape.ColumnCount];
            for (int i = 0; i < def_shape.RowCount; i++)
            {
                for (int j = 0; j < def_shape.ColumnCount; j++)
                {
                    def[i* def_shape.ColumnCount + j] = def_shape[i, j];
                }
            }

            double[] strain = new double[glob_strain.RowCount * glob_strain.ColumnCount];
            double[] stress = new double[glob_stress.RowCount * glob_stress.ColumnCount];
            for (int i = 0; i < glob_stress.RowCount; i++)
            {
                for (int j = 0; j < glob_stress.ColumnCount; j++)
                {
                    stress[i * glob_stress.ColumnCount + j] = glob_stress[i, j];
                    strain[i * glob_stress.ColumnCount + j] = glob_strain[i, j];
                }
            }
            #endregion

            //DA.SetDataTree(0, def_shape_nested);
            //DA.SetDataTree(3, stresses_nested);
            //DA.SetDataTree(4, strain_nested);

            DA.SetDataList(0, def);
            DA.SetDataList(1, reactions);
            DA.SetDataList(2, load);
            DA.SetDataList(3, stress);
            DA.SetDataList(4, strain);
            DA.SetData(5, def_shape);
            DA.SetDataList(6, oldXYZ);

        } //End of main component

        private Grasshopper.DataTree<double> ConvertToNestedList(Matrix<double> M)
        {
            //Convert matrix rows to paths
            Grasshopper.DataTree<double> def_shape_nested = new Grasshopper.DataTree<double>();
            for (int i = 0; i < M.RowCount; i++)
            {
                Grasshopper.Kernel.Data.GH_Path pth = new Grasshopper.Kernel.Data.GH_Path(i);
                for (int j = 0; j < M.ColumnCount; j++)
                {
                    //Adds number to end of current path (pth)
                    def_shape_nested.Add(M[i, j], pth);
                }
            }
            return def_shape_nested;
        }

        private void InterpolateDeformations(Vector<double> def, List<Point3d> points, List<Line> geometry, int n, double height, double width, double alpha, out Matrix<double> def_shape, out List<Point3d> oldXYZ, out Matrix<double> glob_strain)
        {
            def_shape = Matrix<double>.Build.Dense(geometry.Count, (n + 1) * 6);
            glob_strain = Matrix<double>.Build.Dense(geometry.Count, (n + 1) * 3);
            Matrix<double> N, dN;
            Vector<double> u = Vector<double>.Build.Dense(12);
            oldXYZ = new List<Point3d>();
            for (int i = 0; i < geometry.Count; i++)
            {
                //fetches index of original start and endpoint
                Point3d p1 = new Point3d(Math.Round(geometry[i].From.X, 4), Math.Round(geometry[i].From.Y, 4), Math.Round(geometry[i].From.Z, 4));
                Point3d p2 = new Point3d(Math.Round(geometry[i].To.X, 4), Math.Round(geometry[i].To.Y, 4), Math.Round(geometry[i].To.Z, 4));
                int i1 = points.IndexOf(p1);
                int i2 = points.IndexOf(p2);
                //create 12x1 deformation vector for element (6dofs), scaled and populated with existing deformations
                for (int j = 0; j < 6; j++)
                {
                    u[j] = def[i1 * 6 + j];
                    u[j + 6] = def[i2 * 6 + j];
                }

                //interpolate points between startNode and endNode of undeformed (main) element
                List<Point3d> tempOld = InterpolatePoints(geometry[i], n);

                double L = points[i1].DistanceTo(points[i2]);   //L is distance from startnode to endnode

                //Calculate 6 dofs for all new elements using shape functions (n+1 elements)
                Matrix<double> disp = Matrix<double>.Build.Dense(n + 1, 4);
                Matrix<double> rot = Matrix<double>.Build.Dense(n + 1, 4);

                //to show scaled deformations
                Matrix<double> scaled_disp = Matrix<double>.Build.Dense(n + 1, 3);
                
                //transform to local coords
                var tf = TransformationMatrix(geometry[i].From, geometry[i].To, alpha);
                var T = tf.DiagonalStack(tf);
                T = T.DiagonalStack(T);
                u = T * u;

                double x = 0;
                for (int j = 0; j < n + 1; j++)
                {
                    DisplacementField_NB(L, x, out N, out dN);

                    disp.SetRow(j, N.Multiply(u));
                    rot.SetRow(j, dN.Multiply(u));

                    var d0 = new double[] { disp[j, 0], disp[j, 1], disp[j, 2] };
                    var r0 = new double[] { disp[j, 3], rot[j, 2], rot[j, 1] };
                    var t0 = ToGlobal(d0, r0, tf);

                    disp.SetRow(j, new double[] { t0[0], t0[1], t0[2], t0[3] });
                    rot.SetRow(j, new double[] { rot[j, 0], t0[5], t0[4], rot[j, 3] });
                    x += L / n;
                }
                oldXYZ.AddRange(tempOld);

                //add deformation to def_shape (convert from i = nodal number to i = element number)
                def_shape.SetRow(i, SetDef(n + 1, disp, rot));

                glob_strain.SetRow(i, CalculateStrain(n, height, width, u, tf, L, def_shape)); //set strains for all subelement in current element to row i
            }
        }

        private Vector<double> CalculateStrain(int n, double height, double width, Vector<double> u, Matrix<double> tf, double L, Matrix<double> def)
        {
            Matrix<double> dN, ddN;
            double x = 0;
            var strains = Vector<double>.Build.Dense((n + 1) * 3); //contains all subelement strains (only for one element)
            for (int j = 0; j < n + 1; j++)
            {
                DisplacementField_ddN(L, x, out ddN);
                DisplacementField_dN(L, x, out dN);
                
                //u and N are in local coordinates
                var tmp1 = dN * u; //tmp1 = du_x, du_y, du_z, dtheta_x
                var tmp2 = ddN * u; //tmp2 = ddu_x/dx, ddu_y/dx, ddu_z_dx, ddtheta_x/dx
                
                strains[j * 3] = tmp1[0];
                strains[j * 3 + 1] = height * tmp2[2];
                strains[j * 3 + 2] = width * tmp2[1];
                
                x += L / n;
            }
            return strains;
        }

        private Vector<double> ToGlobal(double[] d, double[] r, Matrix<double> tf)
        {
            var dr = Vector<double>.Build.Dense(6);
            for (int i = 0; i < 3; i++)
            {
                dr[i] = d[i];
                dr[i + 3] = r[i];
            }
            tf = tf.DiagonalStack(tf);

            dr = tf.Transpose() * dr;
            return dr;
        }
       
        private double[] SetDef(int m, Matrix<double> disp, Matrix<double> rot)
        {
            //m == n+1
            double[] def_e = new double[m * 6];
            for (int i = 0; i < m; i++)
            {
                //add displacements in x,y,z       
                def_e[i * 6 + 0] = disp[i, 0];
                def_e[i * 6 + 1] = disp[i, 1];
                def_e[i * 6 + 2] = disp[i, 2];
                //add rotations
                def_e[i * 6 + 3] = disp[i, 3];
                def_e[i * 6 + 4] = rot[i, 2]; //theta_y = d_uz/d_x
                def_e[i * 6 + 5] = rot[i, 1]; //theta_z = d_uy/d_x
            }
            return def_e;
        }
        
        private List<Point3d> InterpolatePoints(Line line, int n)
        {
            List<Point3d> tempP = new List<Point3d>(n + 1);
            double[] t = LinSpace(0, 1, n + 1);
            for (int i = 0; i < t.Length; i++)
            {
                var tPm = new Point3d();
                tPm.Interpolate(line.From, line.To, t[i]);
                tPm = new Point3d(Math.Round(tPm.X, 4), Math.Round(tPm.Y, 4), Math.Round(tPm.Z, 4));
                tempP.Add(tPm);
            }
            return tempP;
        }

        private static double[] LinSpace(double x1, double x2, int n)
        {
            //Generate a 1-D array of linearly spaced values
            double step = (x2 - x1) / (n - 1);
            double[] y = new double[n];
            for (int i = 0; i < n; i++)
            {
                y[i] = x1 + step * i;
            }
            return y;
        }

        private void DisplacementField_NB(double L, double x, out Matrix<double> N, out Matrix<double> dN)
        {
            double N1 = 1 - x / L;
            double N2 = x / L;
            double N3 = 1 - 3 * Math.Pow(x, 2) / Math.Pow(L, 2) + 2 * Math.Pow(x, 3) / Math.Pow(L, 3);
            double N4 = x - 2 * Math.Pow(x, 2) / L + Math.Pow(x, 3) / Math.Pow(L, 2);
            double N5 = -N3 + 1;//3 * Math.Pow(x, 2) / Math.Pow(L, 2) - 2 * Math.Pow(x, 3) / Math.Pow(L, 3);
            double N6 = Math.Pow(x, 3) / Math.Pow(L, 2) - Math.Pow(x, 2) / L;

            N = Matrix<double>.Build.DenseOfArray(new double[,] {
                { N1, 0, 0,  0,  0,  0, N2, 0,  0,  0,  0,  0},
                { 0, N3, 0,  0,  0, N4, 0, N5, 0,  0,  0, N6 },
                { 0, 0, N3, 0, -N4, 0, 0, 0, N5, 0, -N6, 0},
                { 0, 0, 0, N1, 0, 0, 0, 0, 0, N2, 0, 0} });

            //u = [u1, u2, u3, u4, u5, u6, u7, u8, u9, u10, u11, u12]
            //u = [ux, uy, uz, theta_x]

            double dN1 = -1 / L;
            double dN2 = 1 / L;
            double dN3 = -6 * x / Math.Pow(L, 2) + 6 * Math.Pow(x, 2) / Math.Pow(L, 3);
            double dN4 = 3 * Math.Pow(x, 2) / Math.Pow(L, 2) - 4 * x / L + 1;
            double dN5 = -dN3;//6 * x / Math.Pow(L, 2) - 6 * Math.Pow(x, 2) / Math.Pow(L, 3);
            double dN6 = 3 * Math.Pow(x, 2) / Math.Pow(L, 2) - 2 * x / L;

            dN = Matrix<double>.Build.DenseOfArray(new double[,] {
                { dN1, 0,    0,  0,  0,     0,  dN2,    0,  0,  0,    0,    0},
                { 0, dN3,    0,  0,  0,  dN4,   0,  dN5,    0,  0,    0, dN6 },
                { 0,    0, dN3,  0, dN4,   0,  0,      0,dN5,  0, dN6,    0},
                //{ 0,    0, dN3,  0, -dN4,   0,  0,      0,dN5,  0, -dN6,    0},
                { 0,    0,   0, dN1, 0,     0,  0,      0,  0, dN2,  0,     0} });

            //theta_y = du_z/dx
            //theta_z = du_y/dx
        }

        private void DisplacementField_dN(double L, double x, out Matrix<double> dN)
        {
            double dN1 = -1 / L;
            double dN2 = 1 / L;
            double dN3 = -6 * x / Math.Pow(L, 2) + 6 * Math.Pow(x, 2) / Math.Pow(L, 3);
            double dN4 = 3 * Math.Pow(x, 2) / Math.Pow(L, 2) - 4 * x / L + 1;
            double dN5 = -dN3;//6 * x / Math.Pow(L, 2) - 6 * Math.Pow(x, 2) / Math.Pow(L, 3);
            double dN6 = 3 * Math.Pow(x, 2) / Math.Pow(L, 2) - 2 * x / L;

            dN = Matrix<double>.Build.DenseOfArray(new double[,] {
                { dN1, 0,    0,  0,  0,     0,  dN2,    0,  0,  0,    0,    0},
                { 0, dN3,    0,  0,  0,  dN4,   0,  dN5,    0,  0,    0, dN6 },
                { 0,    0, dN3,  0, dN4,   0,  0,      0,dN5,  0, dN6,    0},
                //{ 0,    0, dN3,  0, -dN4,   0,  0,      0,dN5,  0, -dN6,    0},
                { 0,    0,   0, dN1, 0,     0,  0,      0,  0, dN2,  0,     0} });

            //theta_y = du_z/dx
            //theta_z = du_y/dx
        }

        private void DisplacementField_ddN(double L, double x, out Matrix<double> ddN)
        {
            double ddN1 = 0;
            double ddN2 = 0;
            double ddN3 = -6 / Math.Pow(L, 2) + 12 * x / Math.Pow(L, 3);
            double ddN4 = -4 / L + 6 * x / Math.Pow(L, 2);
            double ddN5 = 6 / Math.Pow(L, 2) - 12 * x / Math.Pow(L, 3);
            double ddN6 = 6 * x / Math.Pow(L, 2) - 2 / L;

            ddN = Matrix<double>.Build.DenseOfArray(new double[,] {
                { ddN1, 0, 0,  0,  0,  0, ddN2, 0,  0,  0,  0,  0},
                { 0, ddN3, 0,  0,  0, ddN4, 0, ddN5, 0,  0,  0, ddN6 },
                { 0, 0, ddN3, 0, -ddN4,  0, 0, 0, ddN5, 0, -ddN6, 0},
                { 0, 0, 0,  ddN1,  0,  0, 0, 0,  0, ddN2, 0, 0}
            });
        }

        private Vector<double> RestoreTotalDeformationVector(Vector<double> deformations_red, Vector<double> bdc_value)
        {
            Vector<double> def = Vector<double>.Build.Dense(bdc_value.Count);
            for (int i = 0, j = 0; i < bdc_value.Count; i++)
            {
                //if deformation has been calculated, it is added to the vector. Otherwise, the deformation is zero.
                if (bdc_value[i] == 1)
                {
                    def[i] = deformations_red[j];
                    j++;
                }
            }
            return def;
        }

        private void ReducedGlobalStiffnessMatrix(Vector<double> bdc_value, Matrix<double> K, Vector<double> load, out Matrix<double> KGr, out Vector<double> load_red)
        {
            int oldRC = load.Count;
            int newRC = Convert.ToInt16(bdc_value.Sum());
            KGr = Matrix<double>.Build.Dense(newRC, newRC);
            load_red = Vector<double>.Build.Dense(newRC, 0);
            for (int i = 0, ii = 0; i < oldRC; i++)
            {
                //is bdc_value in row i free?
                if (bdc_value[i] == 1)
                {
                    for (int j = 0, jj = 0; j <= i; j++)
                    {
                        //is bdc_value in col j free?
                        if (bdc_value[j] == 1)
                        {
                            //if yes, then add to new K
                            KGr[i - ii, j - jj] = K[i, j];
                            KGr[j - jj, i - ii] = K[i, j];
                        }
                        else
                        {
                            //if not, remember to skip 1 column when adding next time (default matrix value is 0)
                            jj++;
                        }
                    }
                    //add load to reduced list
                    load_red[i - ii] = load[i];
                }
                else
                {
                    //if not, remember to skip 1 row when adding next time (default matrix value is 0)
                    ii++;
                }
            }
        }

        private static bool IsDiagonalPositive(Matrix<double> A)
        {
            bool isPositive = true;
            for (int i = 0; i < A.RowCount; i++)
            {
                if (A[i, i] <= 0)
                {
                    return isPositive = false;
                }
            }
            return isPositive;
        }

        private void CheckSolvers(Matrix<double> K_red, Vector<double> load_red, int decimals, out string output_time)
        {
            string time;
            //Unrounded
            output_time = "Unrounded K_red" + Environment.NewLine;
            bool issym = K_red.IsSymmetric();
            bool diagposi = IsDiagonalPositive(K_red);
            bool ishermitian = K_red.IsHermitian();
            output_time += "IsSymmetric? " + issym + ", IsDiagonalPositive? " + diagposi + ", IsHermitian? " + ishermitian + Environment.NewLine + Environment.NewLine;

            TrySolve(K_red, load_red, out time);
            output_time += "Dense, unrounded K_red: " + Environment.NewLine + time + Environment.NewLine;

            K_red = Matrix<double>.Build.SparseOfMatrix(K_red);
            TrySolve(K_red, load_red, out time);
            output_time += "Sparse, unrounded K_red: " + Environment.NewLine + time + Environment.NewLine;

            //Rounded
            K_red = Matrix<double>.Build.DenseOfMatrix(K_red);
            K_red = RoundMatrix(K_red, decimals);
            output_time += "Rounded to " + decimals + " decimals." + Environment.NewLine;

            issym = K_red.IsSymmetric();
            diagposi = IsDiagonalPositive(K_red);
            ishermitian = K_red.IsHermitian();
            output_time += "IsSymmetric? " + issym + ", IsDiagonalPositive? " + diagposi + ", IsHermitian? " + ishermitian + Environment.NewLine + Environment.NewLine;

            TrySolve(K_red, load_red, out time);
            output_time += "Dense, rounded K_red: " + Environment.NewLine + time + Environment.NewLine;

            K_red = Matrix<double>.Build.SparseOfMatrix(K_red);
            TrySolve(K_red, load_red, out time);
            output_time += "Sparse, rounded K_red: " + Environment.NewLine + time + Environment.NewLine;
            output_time += "=================END OF TEST=================" + Environment.NewLine + Environment.NewLine;
        }

        private void TrySolve(Matrix<double> A, Vector<double> load_red, out string time)
        {
            time = "TrySolve Start:" + Environment.NewLine;
            long timer = 0;
            Stopwatch watch = new Stopwatch();
            try
            {
                watch.Start();
                A.Solve(load_red);
                watch.Stop();
                timer = watch.ElapsedMilliseconds;
                time += "Regular solve: " + timer.ToString() + Environment.NewLine;
            }
            catch (Exception)
            {
                time += "Regular solve raised exception" + Environment.NewLine;
            }

            if (A.GetType().Name == "SparseMatrix")
            {
                time += "End of TrySolve since Matrix is Sparse (other solvers are unsupported)" + Environment.NewLine;
                return;
            }

            try
            {
                watch.Start();
                A.Cholesky().Solve(load_red);
                watch.Stop();
                timer = watch.ElapsedMilliseconds - timer;
                time += "Cholesky solve: " + timer.ToString() + Environment.NewLine;
            }
            catch (Exception)
            {
                time += "Cholesky solve raised exception" + Environment.NewLine;
            }
            try
            {
                watch.Start();
                A.QR().Solve(load_red);
                watch.Stop();
                timer = watch.ElapsedMilliseconds - timer;
                time += "QR solve: " + timer.ToString() + Environment.NewLine;
            }
            catch (Exception)
            {
                time += "QR solve raised exception" + Environment.NewLine;
            }

            try
            {
                watch.Start();
                A.Svd().Solve(load_red);
                watch.Stop();
                timer = watch.ElapsedMilliseconds - timer;
                time += "Svd solve: " + timer.ToString() + Environment.NewLine;
            }
            catch (Exception)
            {
                time += "Svd solve raised exception" + Environment.NewLine;
            }
            try
            {
                watch.Start();
                A.LU().Solve(load_red);
                watch.Stop();
                timer = watch.ElapsedMilliseconds - timer;
                time += "LU solve: " + timer.ToString() + Environment.NewLine;
            }
            catch (Exception)
            {
                time += "LU solve raised exception" + Environment.NewLine;
            }
        }

        private string PrintMatrix(object B, string head)
        {
            Debug.WriteLine(head);
            string ss = "";
            if (B.GetType() == typeof(SparseMatrix))
            {
                Matrix<double> A = (Matrix<double>)B;
                ss += A.RowCount.ToString() + "x" + A.ColumnCount.ToString() + "-matrix" + Environment.NewLine;
                for (int i = 0; i < A.RowCount; i++)
                {
                    for (int j = 0; j < A.ColumnCount; j++)
                    {
                        ss += String.Format(" {0,15:0.0} ", A[i, j]);
                    }
                    ss += Environment.NewLine + Environment.NewLine + Environment.NewLine;
                }
                Debug.WriteLine(ss);
                return ss;
            }
            else if (B.GetType() == typeof(DenseMatrix))
            {
                Matrix<double> A = (Matrix<double>)B;
                ss += A.RowCount.ToString() + "x" + A.ColumnCount.ToString() + "-matrix" + Environment.NewLine;
                for (int i = 0; i < A.RowCount; i++)
                {
                    for (int j = 0; j < A.ColumnCount; j++)
                    {
                        ss += String.Format(" {0,15:0.0} ", A[i, j]);
                    }
                    ss += Environment.NewLine + Environment.NewLine + Environment.NewLine;
                }
                Debug.WriteLine(ss);
                return ss;
            }
            else if (B.GetType() == typeof(MathNet.Numerics.LinearAlgebra.Double.DenseVector))
            {
                Vector<double> A = (Vector<double>)B;
                ss += "1x" + A.Count.ToString() + "-vector" + Environment.NewLine;
                for (int i = 0; i < A.Count; i++)
                {
                    ss += String.Format(" {0,10:0.0} ", A[i]);
                }
                System.Diagnostics.Debug.WriteLine(ss);
                return ss;
            }
            else if (B.GetType() == typeof(List<int>))
            {
                List<int> A = (List<int>)B;
                ss += "1x" + A.Count.ToString() + "-list" + Environment.NewLine;
                for (int i = 0; i < A.Count; i++)
                {
                    ss += String.Format(" {0,10:0.0} ", A[i]);
                }
                System.Diagnostics.Debug.WriteLine(ss);
                return ss;
            }
            else if (B.GetType() == typeof(List<double>))
            {
                List<double> A = (List<double>)B;
                ss += "1x" + A.Count.ToString() + "-list" + Environment.NewLine;
                for (int i = 0; i < A.Count; i++)
                {
                    ss += String.Format(" {0,10:0.0} ", A[i]);
                }
                System.Diagnostics.Debug.WriteLine(ss);
                return ss;
            }
            return ss;
        }

        private static Matrix<double> RoundMatrix(Matrix<double> A, int decs)
        {
            for (int i = 0; i < A.RowCount; i++)
            {
                for (int j = 0; j < A.ColumnCount; j++)
                {
                    A[i, j] = Math.Round(A[i, j], decs);
                }
            }
            return A;
        }

        private Matrix<double> TransformationMatrix(Point3d p1, Point3d p2, double alpha)
        {
            double L = p1.DistanceTo(p2);

            double cx = (p2.X - p1.X) / L;
            double cy = (p2.Y - p1.Y) / L;
            double cz = (p2.Z - p1.Z) / L;
            double c1 = Math.Cos(alpha);
            double s1 = Math.Sin(alpha);
            double cxz = Math.Round(Math.Sqrt(Math.Pow(cx, 2) + Math.Pow(cz, 2)), 6);

            Matrix<double> t;

            if (Math.Round(cx, 6) == 0 && Math.Round(cz, 6) == 0)
            {
                t = Matrix<double>.Build.DenseOfArray(new double[,]
            {
                    {      0, cy,  0},
                    { -cy*c1,  0, s1},
                    {  cy*s1,  0, c1},
            });
            }
            else
            {
                t = Matrix<double>.Build.DenseOfArray(new double[,]
            {
                    {                     cx,       cy,                   cz},
                    {(-cx*cy*c1 - cz*s1)/cxz,   cxz*c1,(-cy*cz*c1+cx*s1)/cxz},
                    {   (cx*cy*s1-cz*c1)/cxz,  -cxz*s1, (cy*cz*s1+cx*c1)/cxz},
            });
            }
            return t;
        }

        private void ElementStiffnessMatrix(Line currentLine, double E, double A, double Iy, double Iz, double J, double G, double alpha, out Point3d p1, out Point3d p2, out Matrix<double> Ke)
        {
            double L = Math.Round(currentLine.Length, 6);

            p1 = new Point3d(Math.Round(currentLine.From.X, 4), Math.Round(currentLine.From.Y, 4), Math.Round(currentLine.From.Z, 4));
            p2 = new Point3d(Math.Round(currentLine.To.X, 4), Math.Round(currentLine.To.Y, 4), Math.Round(currentLine.To.Z, 4));

            Matrix<double> tf = TransformationMatrix(p1, p2, alpha);
            var T = tf.DiagonalStack(tf);
            T = T.DiagonalStack(T);

            Matrix<double> T_T = T.Transpose();

            double A1 = (E * A) / (L);

            double kz1 = (12 * E * Iz) / (L * L * L);
            double kz2 = (6 * E * Iz) / (L * L);
            double kz3 = (4 * E * Iz) / L;
            double kz4 = (2 * E * Iz) / L;

            double ky1 = (12 * E * Iy) / (L * L * L);
            double ky2 = (6 * E * Iy) / (L * L);
            double ky3 = (4 * E * Iy) / L;
            double ky4 = (2 * E * Iy) / L;

            double C1 = (G * J) / L;

            Matrix<double> ke = DenseMatrix.OfArray(new double[,]
            {
                    { A1,    0,    0,    0,    0,    0,  -A1,    0,    0,    0,    0,    0 },
                    {  0,  kz1,    0,    0,    0,  kz2,    0, -kz1,    0,    0,    0,  kz2 },
                    {  0,    0,  ky1,    0, -ky2,    0,    0,    0, -ky1,    0, -ky2,    0 },
                    {  0,    0,    0,   C1,    0,    0,    0,    0,    0,  -C1,    0,    0 },
                    {  0,    0, -ky2,    0,  ky3,    0,    0,    0,  ky2,    0,  ky4,    0 },
                    {  0,  kz2,    0,    0,    0,  kz3,    0, -kz2,    0,    0,    0,  kz4 },
                    {-A1,    0,    0,    0,    0,    0,   A1,    0,    0,    0,    0,    0 },
                    {  0, -kz1,    0,    0,    0, -kz2,    0,  kz1,    0,    0,    0, -kz2 },
                    {  0,    0, -ky1,    0,  ky2,    0,    0,    0,  ky1,    0,  ky2,    0 },
                    {  0,    0,    0,  -C1,    0,    0,    0,    0,    0,   C1,    0,    0 },
                    {  0,    0, -ky2,    0,  ky4,    0,    0,    0,  ky2,    0,  ky3,    0 },
                    {  0,  kz2,    0,    0,    0,  kz4,    0, -kz2,    0,    0,    0,  kz3 },
            });

            ke = ke.Multiply(T);
            Ke = T_T.Multiply(ke);
        }

        private Matrix<double> GlobalStiffnessMatrix(List<Line> geometry, List<Point3d> points, double E, double A, double Iy, double Iz, double J, double G, double alpha)
        {
            int gdofs = points.Count * 6;
            Matrix<double> KG = DenseMatrix.OfArray(new double[gdofs, gdofs]);

            foreach (Line currentLine in geometry)
            {
                Matrix<double> Ke;
                Point3d p1, p2;

                //Calculate Ke
                ElementStiffnessMatrix(currentLine, E, A, Iy, Iz, J, G, alpha, out p1, out p2, out Ke);

                //Fetch correct point indices
                int node1 = points.IndexOf(p1);
                int node2 = points.IndexOf(p2);

                //Inputting Ke to correct entries in Global Stiffness Matrix
                for (int i = 0; i < Ke.RowCount / 2; i++)
                {
                    for (int j = 0; j < Ke.ColumnCount / 2; j++)
                    {
                        //top left 3x3 of k-element matrix
                        KG[node1 * 6 + i, node1 * 6 + j] += Ke[i, j];
                        //top right 3x3 of k-element matrix  
                        KG[node1 * 6 + i, node2 * 6 + j] += Ke[i, j + 6];
                        //bottom left 3x3 of k-element matrix
                        KG[node2 * 6 + i, node1 * 6 + j] += Ke[i + 6, j];
                        //bottom right 3x3 of k-element matrix
                        KG[node2 * 6 + i, node2 * 6 + j] += Ke[i + 6, j + 6];
                    }
                }
            }
            return KG;
        }

        private Vector<double> CreateLoadList(List<string> loadtxt, List<string> momenttxt, List<Point3d> points)
        {
            Vector<double> loads = Vector<double>.Build.Dense(points.Count * 6);
            List<double> inputLoads = new List<double>();
            List<Point3d> coordlist = new List<Point3d>();

            for (int i = 0; i < loadtxt.Count; i++)
            {
                string coordstr = (loadtxt[i].Split(':')[0]);
                string loadstr = (loadtxt[i].Split(':')[1]);

                string[] coordstr1 = (coordstr.Split(','));
                string[] loadstr1 = (loadstr.Split(','));

                inputLoads.Add(Math.Round(double.Parse(loadstr1[0]), 4));
                inputLoads.Add(Math.Round(double.Parse(loadstr1[1]), 4));
                inputLoads.Add(Math.Round(double.Parse(loadstr1[2]), 4));

                coordlist.Add(new Point3d(Math.Round(double.Parse(coordstr1[0]), 4), Math.Round(double.Parse(coordstr1[1]), 4), Math.Round(double.Parse(coordstr1[2]), 4)));
            }

            foreach (Point3d point in coordlist)
            {
                int i = points.IndexOf(point);
                int j = coordlist.IndexOf(point);
                loads[i * 6 + 0] = inputLoads[j * 3 + 0]; //is loads out of range? (doesn't seem to have been initialized with size yet)
                loads[i * 6 + 1] = inputLoads[j * 3 + 1];
                loads[i * 6 + 2] = inputLoads[j * 3 + 2];
            }
            inputLoads.Clear();
            coordlist.Clear();
            for (int i = 0; i < momenttxt.Count; i++) if (momenttxt[0] != "")
                {
                    string coordstr = (momenttxt[i].Split(':')[0]);
                    string loadstr = (momenttxt[i].Split(':')[1]);

                    string[] coordstr1 = (coordstr.Split(','));
                    string[] loadstr1 = (loadstr.Split(','));

                    inputLoads.Add(Math.Round(double.Parse(loadstr1[0]), 4));
                    inputLoads.Add(Math.Round(double.Parse(loadstr1[1]), 4));
                    inputLoads.Add(Math.Round(double.Parse(loadstr1[2]), 4));


                    coordlist.Add(new Point3d(Math.Round(double.Parse(coordstr1[0]), 4), Math.Round(double.Parse(coordstr1[1]), 4), Math.Round(double.Parse(coordstr1[2]), 4)));
                }

            foreach (Point3d point in coordlist)
            {
                int i = points.IndexOf(point);
                int j = coordlist.IndexOf(point);
                loads[i * 6 + 3] = inputLoads[j * 3 + 0];
                loads[i * 6 + 4] = inputLoads[j * 3 + 1];
                loads[i * 6 + 5] = inputLoads[j * 3 + 2];
            }
            return loads;
        }

        private Vector<double> CreateBDCList(List<string> bdctxt, List<Point3d> points)
        {
            //initializing bdc_value as vector of size gdofs, and entry values = 1
            Vector<double> bdc_value = Vector.Build.Dense(points.Count * 6, 1);
            List<int> bdcs = new List<int>();
            List<Point3d> bdc_points = new List<Point3d>(); //Coordinates relating til bdc_value in for (eg. x y z)

            //Parse string input
            for (int i = 0; i < bdctxt.Count; i++)
            {
                string coordstr = (bdctxt[i].Split(':')[0]);
                string bdcstr = (bdctxt[i].Split(':')[1]);

                string[] coordstr1 = (coordstr.Split(','));
                string[] bdcstr1 = (bdcstr.Split(','));

                bdc_points.Add(new Point3d(Math.Round(double.Parse(coordstr1[0]), 4), Math.Round(double.Parse(coordstr1[1]), 4), Math.Round(double.Parse(coordstr1[2]), 4)));

                bdcs.Add(int.Parse(bdcstr1[0]));
                bdcs.Add(int.Parse(bdcstr1[1]));
                bdcs.Add(int.Parse(bdcstr1[2]));
                bdcs.Add(int.Parse(bdcstr1[3]));
                bdcs.Add(int.Parse(bdcstr1[4]));
                bdcs.Add(int.Parse(bdcstr1[5]));
            }

            //Format to correct entries in bdc_value
            foreach (var point in bdc_points)
            {
                int globalI = points.IndexOf(point);
                int localI = bdc_points.IndexOf(point);
                bdc_value[globalI * 6 + 0] = bdcs[localI * 6 + 0];
                bdc_value[globalI * 6 + 1] = bdcs[localI * 6 + 1];
                bdc_value[globalI * 6 + 2] = bdcs[localI * 6 + 2];
                bdc_value[globalI * 6 + 3] = bdcs[localI * 6 + 3];
                bdc_value[globalI * 6 + 4] = bdcs[localI * 6 + 4];
                bdc_value[globalI * 6 + 5] = bdcs[localI * 6 + 5];
            }
            return bdc_value;
        }
        
        private void SetMaterial(string mattxt, out double E, out double A, out double Iy, out double Iz, out double J, out double G, out double v, out double alpha)
        {
            string[] matProp = (mattxt.Split(','));

            E = (Math.Round(double.Parse(matProp[0]), 2));
            A = (Math.Round(double.Parse(matProp[1]), 2));
            Iy = (Math.Round(double.Parse(matProp[2]), 2));
            Iz = (Math.Round(double.Parse(matProp[3]), 2));
            v = (Math.Round(double.Parse(matProp[4]), 2));
            G = E / (2 * (1 + Math.Pow(v, 2)));
            alpha = (Math.Round(double.Parse(matProp[5]), 2))*Math.PI/180; //to radians

            J = Iy + Iz;
        }

        private List<Point3d> CreatePointList(List<Line> geometry)
        {
            List<Point3d> points = new List<Point3d>();
            foreach (Line line in geometry) //adds point unless it already exists in pointlist
            {
                Point3d tempFrom = new Point3d(Math.Round(line.From.X, 4), Math.Round(line.From.Y, 4), Math.Round(line.From.Z, 4));
                Point3d tempTo = new Point3d(Math.Round(line.To.X, 4), Math.Round(line.To.Y, 4), Math.Round(line.To.Z, 4));

                if (!points.Contains(tempFrom))
                {
                    points.Add(tempFrom);
                }
                if (!points.Contains(tempTo))
                {
                    points.Add(tempTo);
                }
            }
            return points;
        }

        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                return Properties.Resources.Calc;
            }
        }

        public override Guid ComponentGuid
        {
            get { return new Guid("d636ebc9-0d19-44d5-a3ad-cec704b82323"); }
        }


        /// Component Visual//
        public class Attributes_Custom : Grasshopper.Kernel.Attributes.GH_ComponentAttributes
        {
            public Attributes_Custom(GH_Component owner) : base(owner) { }
            protected override void Layout()
            {
                base.Layout();

                Rectangle rec0 = GH_Convert.ToRectangle(Bounds);

                rec0.Height += 22;

                Rectangle rec1 = rec0;
                rec1.X = rec0.Left + 1;
                rec1.Y = rec0.Bottom - 22;
                rec1.Width = (rec0.Width) / 3 + 1;
                rec1.Height = 22;
                rec1.Inflate(-2, -2);

                //Rectangle rec2 = rec1;
                //rec2.X = rec1.Right + 2;

                //Rectangle rec3 = rec1;
                //rec3.X = rec2.Right + 2;

                Bounds = rec0;
                ButtonBounds = rec1;
                //ButtonBounds2 = rec2;
                //ButtonBounds3 = rec3;

            }

            GH_Palette xColor = GH_Palette.Grey;
            //GH_Palette yColor = GH_Palette.Grey;
            //GH_Palette zColor = GH_Palette.Grey;

            private Rectangle ButtonBounds { get; set; }
            //private Rectangle ButtonBounds2 { get; set; }
            //private Rectangle ButtonBounds3 { get; set; }

            protected override void Render(GH_Canvas canvas, Graphics graphics, GH_CanvasChannel channel)
            {
                base.Render(canvas, graphics, channel);
                if (channel == GH_CanvasChannel.Objects)
                {
                    GH_Capsule button = GH_Capsule.CreateTextCapsule(ButtonBounds, ButtonBounds, xColor, "Run", 3, 0);
                    button.Render(graphics, Selected, false, false);
                    button.Dispose();
                }
                //if (channel == GH_CanvasChannel.Objects)
                //{
                //    GH_Capsule button2 = GH_Capsule.CreateTextCapsule(ButtonBounds2, ButtonBounds2, yColor, "Run Test", 2, 0);
                //    button2.Render(graphics, Selected, Owner.Locked, false);
                //    button2.Dispose();
                //}
                //if (channel == GH_CanvasChannel.Objects)
                //{
                //    GH_Capsule button3 = GH_Capsule.CreateTextCapsule(ButtonBounds3, ButtonBounds3, zColor, "AboveLimit", 2, 0);
                //    button3.Render(graphics, Selected, Owner.Locked, false);
                //    button3.Dispose();
                //}
            }

            public override GH_ObjectResponse RespondToMouseDown(GH_Canvas sender, GH_CanvasMouseEvent e)
            {
                if (e.Button == MouseButtons.Left)
                {
                    RectangleF rec = ButtonBounds;
                    if (rec.Contains(e.CanvasLocation))
                    {
                        switchColor("Run");
                        if (xColor == GH_Palette.Black) { CalcComponent.setStart("Run", true); Owner.ExpireSolution(true); }
                        if (xColor == GH_Palette.Grey) { CalcComponent.setStart("Run", false); }
                        sender.Refresh();
                        return GH_ObjectResponse.Handled;
                    }
                    //rec = ButtonBounds2;
                    //if (rec.Contains(e.CanvasLocation))
                    //{
                    //    switchColor("Run Test");
                    //    if (yColor == GH_Palette.Black) { CalcComponent.setStart("Run Test", true); Owner.ExpireSolution(true); }
                    //    if (yColor == GH_Palette.Grey) { CalcComponent.setStart("Run Test", false); }
                    //    sender.Refresh();
                    //    return GH_ObjectResponse.Handled;
                    //}
                    //rec = ButtonBounds3;
                    //if (rec.Contains(e.CanvasLocation))
                    //{
                    //    switchColor("AboveLimit");
                    //    if (zColor == GH_Palette.Black) { CalcComponent.setStart("AboveLimit", true); Owner.ExpireSolution(true); }
                    //    if (zColor == GH_Palette.Grey) { CalcComponent.setStart("AboveLimit", false); }
                    //    sender.Refresh();
                    //    return GH_ObjectResponse.Handled;
                    //}
                }
                return base.RespondToMouseDown(sender, e);
            }

            private void switchColor(string button)
            {
                if (button == "Run")
                {
                    if (xColor == GH_Palette.Black) { xColor = GH_Palette.Grey; }
                    else { xColor = GH_Palette.Black; }
                }
                //else if (button == "Run Test")
                //{
                //    if (yColor == GH_Palette.Black) { yColor = GH_Palette.Grey; }
                //    else { yColor = GH_Palette.Black; }
                //}
                //else if (button == "AboveLimit")
                //{
                //    if (zColor == GH_Palette.Black) { zColor = GH_Palette.Grey; }
                //    else { zColor = GH_Palette.Black; }
                //}
            }
        }
    }
}

//private void CalculateStrains(List<Point3d> old, List<Point3d> new1, int n, int m, double z, Matrix<double> def_shape, out Matrix<double> glob_strain)
//{
//    glob_strain = Matrix<double>.Build.Dense(m, n); //geometry.Count rows, n columns
//    for (int i = 0; i < m; i++)
//    {
//        for (int j = 0; j < n; j++)
//        {
//            var originL = old[m + j].DistanceTo(old[m + j + 1]);
//            var dL = new1[m + j].DistanceTo(new1[m + j + 1]) - originL;
//            var eps_x_ax = dL / originL;

//            TransformationMatrix()
//            var eps_x_b = 
//        }
//    }
//}

//private List<string> AboveStressLimit(Matrix<double> mises_stress, double limit)
//{
//    List<string> s = new List<string>();
//    for (int i = 0; i < mises_stress.RowCount; i++)
//    {
//        string ts = "";
//        string ds = "";
//        bool f = true;
//        for (int j = 0; j < mises_stress.ColumnCount; j++)
//        {
//            if (mises_stress[i, j] > limit)
//            {
//                if (f)
//                {
//                    ts += j;
//                    ds += Math.Round(mises_stress[i, j]);
//                    f = false;
//                }
//                else
//                {
//                    ts += ", " + j;
//                    ds += ", " + Math.Round(mises_stress[i, j]);
//                }
//            }
//        }
//        if (ts.Length > 0)
//        {
//            s.Add("Main element: " + i + "; Sub element(s): " + ts + "; Stresses: " + ds);
//        }
//    }
//    if (s.Count == 0)
//    {
//        s.Add("No elements had stresses above: " + limit + "[MPa]");
//    }
//    return s;
//}

//private void CalculateMisesStresses(Matrix<double> glob_stress, out Matrix<double> mises_stress)
//{
//    double sig_v;
//    mises_stress = Matrix<double>.Build.Dense(glob_stress.RowCount, glob_stress.ColumnCount / 6);
//    for (int i = 0; i < glob_stress.RowCount; i++)
//    {
//        for (int j = 0; j < glob_stress.ColumnCount/6; j++)
//        {
//            sig_v = Math.Sqrt(0.5 *(
//                Math.Pow(glob_stress[i, j * 6] - glob_stress[i, j * 6 + 1], 2) + 
//                Math.Pow(glob_stress[i, j * 6 + 1] - glob_stress[i, j * 6 + 2], 2) + 
//                Math.Pow(glob_stress[i, j * 6 + 2] - glob_stress[i, j * 6], 2)
//               + 3 * ( //NB! multiplication factor 3 is for shear stress12, 23 and 31, while mf should be 6 for stress12, 23, 31
//               glob_stress[i, j * 6 + 3] * glob_stress[i, j * 6 + 3] + 
//               glob_stress[i, j * 6 + 4] * glob_stress[i, j * 6 + 4] + 
//               glob_stress[i, j * 6 + 5] * glob_stress[i, j * 6 + 5]))); 
//            mises_stress[i, j] = sig_v;
//        }
//    }
//}
//
//private void CalculateInternalStresses(Matrix<double> strain, double E, double G, double v, out Matrix<double> stress)
//{                
//    //strain = [ex, ey, ez, gxy, gyz, gzx]
//    stress = Matrix<double>.Build.Dense(strain.RowCount, strain.ColumnCount);
//    double c = E / ((1.0 + v)*(1.0 - 2.0*v));
//    for (int i = 0; i < strain.RowCount; i++)
//    {
//        Debug.WriteLine(strain.Row(i));
//        for (int j = 0; j < strain.ColumnCount/6; j++)
//        {
//            stress[i, j * 6 + 0] = c * ((1 - v) * strain[i, j * 6 + 0] + v * (strain[i, j * 6 + 1] + strain[i, j * 6 + 2]));
//            stress[i, j * 6 + 1] = c * ((1 - v) * strain[i, j * 6 + 1] + v * (strain[i, j * 6 + 0] + strain[i, j * 6 + 2]));
//            stress[i, j * 6 + 2] = c * ((1 - v) * strain[i, j * 6 + 2] + v * (strain[i, j * 6 + 0] + strain[i, j * 6 + 1]));
//            stress[i, j * 6 + 3] = G * strain[i, j * 6 + 3];
//            stress[i, j * 6 + 4] = G * strain[i, j * 6 + 4];
//            stress[i, j * 6 + 5] = G * strain[i, j * 6 + 5];
//            Debug.WriteLine(stress.Row(i));
//        }
//    }
//}

#region Old CalcComp
//using System;
//using System.Collections.Generic;

//using Grasshopper.Kernel;
//using Rhino.Geometry;
//using System.Drawing;
//using Grasshopper.GUI.Canvas;
//using System.Windows.Forms;
//using Grasshopper.GUI;

//using MathNet.Numerics.LinearAlgebra;
//using MathNet.Numerics.LinearAlgebra.Double;
//using System.Diagnostics;

//namespace Beam3D
//{
//    public class CalcComponent : GH_Component
//    {
//        public CalcComponent()
//          : base("BeamCalculation", "BeamC",
//              "Description",
//              "Koala", "3D Beam")
//        {
//        }

//        //Initialize moments
//        static bool startCalc = true;
//        static bool startTest = false;

//        //Method to allow c hanging of variables via GUI (see Component Visual)
//        public static void setStart(string s, bool i)
//        {
//            if (s == "Run")
//            {
//                startCalc = i;
//            }
//            if (s == "Run Test")
//            {
//                startTest = i;
//            }
//            Grasshopper.Instances.ActiveCanvas.Document.ExpireSolution();
//            Grasshopper.Instances.ActiveCanvas.Document.NewSolution(false);
//        }

//        public override void CreateAttributes()
//        {
//            m_attributes = new Attributes_Custom(this);
//        }

//        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
//        {
//            pManager.AddLineParameter("Lines", "LNS", "Geometry, in form of Lines)", GH_ParamAccess.list);
//            pManager.AddTextParameter("Boundary Conditions", "BDC", "Boundary Conditions in form x,y,z,vx,vy,vz,rx,ry,rz", GH_ParamAccess.list);
//            pManager.AddTextParameter("Material properties", "Mat", "Material Properties", GH_ParamAccess.item, "210000,3600,4920000,4920000,79300");
//            pManager.AddTextParameter("PointLoads", "PL", "Load given as Vector [N]", GH_ParamAccess.list);
//            pManager.AddTextParameter("PointMoment", "PM", "Moment set in a point in [Nm]", GH_ParamAccess.list, "");
//            pManager.AddIntegerParameter("Elements", "n", "Number of elements", GH_ParamAccess.item, 1);
//        }

//        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
//        {
//            pManager.AddNumberParameter("Deformations", "Def", "Deformations", GH_ParamAccess.list);
//            pManager.AddNumberParameter("Reactions", "R", "Reaction Forces", GH_ParamAccess.list);
//            pManager.AddNumberParameter("Element stresses", "Strs", "The Stress in each element", GH_ParamAccess.list);
//            pManager.AddNumberParameter("Element strains", "Strn", "The Strain in each element", GH_ParamAccess.list);
//        }

//        protected override void SolveInstance(IGH_DataAccess DA)
//        {
//            #region Fetch input
//            //Expected inputs
//            List<Line> geometry = new List<Line>();         //Initial Geometry of lines
//            List<string> bdctxt = new List<string>();       //Boundary conditions in string format
//            List<string> loadtxt = new List<string>();      //loads in string format
//            List<string> momenttxt = new List<string>();    //Moments in string format
//            string mattxt = "";
//            int n = 1;


//            //Set expected inputs from Indata
//            if (!DA.GetDataList(0, geometry)) return;       //sets geometry
//            if (!DA.GetDataList(1, bdctxt)) return;         //sets boundary conditions as string
//            if (!DA.GetData(2, ref mattxt)) return;         //sets material properties as string
//            if (!DA.GetDataList(3, loadtxt)) return;        //sets load as string
//            if (!DA.GetDataList(4, momenttxt)) return;      //sets moment as string
//            if (!DA.GetData(5, ref n)) return;              //sets number of elements
//            #endregion

//            //Interpret and set material parameters
//            double E;       //Material Young's modulus, initial value 210000 [MPa]
//            double A;       //Area for each element in same order as geometry, initial value CFS100x100 3600 [mm^2]
//            double Iy;      //Moment of inertia about local y axis, initial value 4.92E6 [mm^4]
//            double Iz;      //Moment of inertia about local z axis, initial value 4.92E6 [mm^4]
//            double J;       //Polar moment of inertia
//            double G;       //Shear modulus, initial value 79300 [mm^4]
//            SetMaterial(mattxt, out E, out A, out Iy, out Iz, out J, out G);

//            #region Prepares geometry, boundary conditions and loads for calculation
//            //List all nodes (every node only once), numbering them according to list index
//            List<Point3d> points = CreatePointList(geometry);


//            //Interpret the BDC inputs (text) and create list of boundary condition (1/0 = free/clamped) for each dof.
//            List<int> bdc_value = CreateBDCList(bdctxt, points);


//            //Interpreting input load (text) and creating load list (do uble)
//            List<double> load = CreateLoadList(loadtxt, momenttxt, points);
//            #endregion

//            Vector<double> def_tot;
//            Vector<double> reactions;
//            List<double> internalStresses;
//            List<double> internalStrains;

//            if (startCalc)
//            {
//                #region Create global and reduced stiffness matrix
//                //Create global stiffness matrix
//                Matrix<double> K_tot = GlobalStiffnessMatrix(geometry, points, E, A, Iy, Iz, J, G);

//                //Create reduced K-matrix and reduced load list (removed free dofs)
//                Matrix<double> K_red;
//                Vector<double> load_red;
//                CreateReducedGlobalStiffnessMatrix(bdc_value, K_tot, load, out K_red, out load_red);
//                #endregion

//                #region Solver Performance Test
//                if (startTest)
//                {
//                    string output_time = "";
//                    string performanceResult = "=================START OF TEST=================" + Environment.NewLine;
//                    performanceResult += "Number of lines: " + geometry.Count.ToString() + Environment.NewLine;
//                    string tester = "";

//                    //checking for error with writing to file (skip test if unable to write)
//                    try
//                    {
//                        using (System.IO.StreamWriter file =
//                        new System.IO.StreamWriter(@"solverTest.txt", true))
//                        {
//                            file.WriteLine(tester);
//                        }
//                    }
//                    //create new file if solverTest.txt does not exist
//                    catch (System.IO.DirectoryNotFoundException)
//                    {
//                        System.IO.File.WriteAllText(@"solverTest.txt", tester);
//                    }
//                    //other exception (no write access?)
//                    catch (Exception)
//                    {
//                        AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Write to file error! Create file solverTest.txt in Koala folder and/or Rhinoceros 5.exe folder");
//                        return;
//                    }

//                    int decimals = 2;
//                    CheckSolvers(K_red, load_red, decimals, out output_time);

//                    performanceResult += output_time;

//                    //append result to txt-file (at buildpath)
//                    try
//                    {
//                        using (System.IO.StreamWriter file =
//                        new System.IO.StreamWriter(@"solverTest.txt", true))
//                        {
//                            file.WriteLine(performanceResult);
//                        }
//                    }
//                    //create new file if solverTest.txt does not exist
//                    catch (System.IO.DirectoryNotFoundException)
//                    {
//                        System.IO.File.WriteAllText(@"\solverTest.txt", output_time);
//                    }

//                }
//                #endregion

//                #region Calculate deformations, reaction forces and internal strains and stresses
//                //Calculate deformations
//                Vector<double> def_red = K_red.Cholesky().Solve(load_red);

//                //Add the clamped dofs (= 0) to the deformations list
//                def_tot = RestoreTotalDeformationVector(def_red, bdc_value);

//                //Interpolate deformations using shape functions
//                Vector<double> def_shape = InterpolateDeformations(def_tot, points, geometry, n);

//                //Calculate the reaction forces from the deformations
//                reactions = K_tot.Multiply(def_tot);

//                //Calculate the internal strains and stresses in each member
//                CalculateInternalStrainsAndStresses(def_shape, points, E, geometry, out internalStresses, out internalStrains);
//                #endregion
//            }
//            else
//            {
//                def_tot = Vector<double>.Build.Dense(points.Count * 6);
//                reactions = def_tot;

//                internalStresses = new List<double>(geometry.Count);
//                internalStresses.AddRange(new double[geometry.Count]);
//                internalStrains = internalStresses;
//            }
//            DA.SetDataList(0, def_tot);
//            DA.SetDataList(1, reactions);
//            DA.SetDataList(2, internalStresses);
//            DA.SetDataList(3, internalStrains);


//        } //End of main component

//        private Vector<double> InterpolateDeformations(Vector<double> def, List<Point3d> points, List<Line> geometry, int n)
//        {
//            List<Point3d> newXYZ = new List<Point3d>();
//            List<Point3d> oldXYZ = new List<Point3d>();
//            Matrix<double> rotations = Matrix<double>.Build.Dense(2 * n * 3, 4);
//            Vector<double> shape_def = Vector<double>.Build.Dense(2 * n * 6);
//            Matrix<double> N, B;

//            int counter = 0;
//            foreach (Line line in geometry)
//            {
//                //fetches index of original start and endpoint
//                int i1 = points.IndexOf(line.From);
//                int i2 = points.IndexOf(line.To);

//                //create 12x1 deformation vector for element (6dofs), scaled and populated with existing deformations
//                var u = Vector<double>.Build.Dense(12);
//                for (int j = 0; j < 6; j++)
//                {
//                    u[j] = def[i1 * 6 + j];
//                    u[j + 6] = def[i2 * 6 + j];
//                }
//                //u = scale * u;

//                //interpolate points between line.From and line.To
//                List<Point3d> tempP = InterpolatePoints(line, n);
//                oldXYZ.AddRange(tempP);

//                double L = points[i1].DistanceTo(points[i2]);   //L is distance from startnode to endnode
//                var x = Vector<double>.Build.Dense(n + 1);      //maybe this should be projected x instead???
//                for (int j = 0; j < n + 1; j++)
//                {
//                    x[j] = j * L / n;
//                }


//                //Calculate 6 dofs for all new elements using shape functions (n+1 elements)
//                Matrix<double> disp = Matrix<double>.Build.Dense(n + 1, 4);
//                Matrix<double> rot = Matrix<double>.Build.Dense(n + 1, 4);

//                for (int j = 0; j < n + 1; j++)          //x are points inbetween (?)
//                {
//                    Shapefunctions(L, x[j], out N, out B);
//                    disp.SetRow(j, N.Multiply(u));
//                    rot.SetRow(j, B.Multiply(u));
//                    rotations.SetRow(counter + j, rot.Row(j));
//                    counter++;
//                }

//                //Calculate new nodal points
//                for (int j = 0; j < n + 1; j++)
//                {
//                    //original xyz                        
//                    var tP = tempP[j];

//                    //add displacement
//                    //tP.X += disp[j, 0];
//                    //tP.Y += disp[j, 1];
//                    //tP.Z += disp[j, 2];

//                    tP.X = tP.X + disp[j, 0] + tP.Z * Math.Cos(Math.PI / 2 - rot[j, 2]) + tP.Y * Math.Cos(Math.PI / 2 - rot[j, 1]);
//                    tP.Y = -Math.Cos(disp[j, 3]) * tP.Y * Math.Sin(Math.PI / 2 - rot[j, 1]) + Math.Sin(disp[j, 3]) * tP.Z + disp[j, 1];
//                    //tP.Z = disp[j, 2] + tP.Z * Math.Sin(rot[j, 2]);
//                    tP.Z = -Math.Sin(disp[j, 3]) * tP.Y + Math.Cos(disp[j, 3]) * tP.Z * Math.Sin(Math.PI / 2 - rot[j, 3]) + disp[j, 2];

//                    //replace previous xyz with displaced xyz
//                    tempP[j] = tP;
//                }
//                newXYZ.AddRange(tempP);
//            }
//            //calculate distance from original interpolated points to deformed points
//            for (int i = 0; i < oldXYZ.Count; i++)
//            {
//                shape_def[i * 6 + 0] = newXYZ[i].X - oldXYZ[i].X;
//                shape_def[i * 6 + 1] = newXYZ[i].Y - oldXYZ[i].Y;
//                shape_def[i * 6 + 2] = newXYZ[i].Z - oldXYZ[i].Z;
//                shape_def[i * 6 + 3] = rotations[i, 0];
//                shape_def[i * 6 + 4] = rotations[i, 1];
//                shape_def[i * 6 + 5] = rotations[i, 2];
//            }
//            return shape_def;
//        }

//        private List<Point3d> InterpolatePoints(Line line, int n)
//        {
//            List<Point3d> tempP = new List<Point3d>(n + 1);
//            double[] t = LinSpace(0, 1, n + 1);
//            for (int i = 0; i < t.Length; i++)
//            {
//                var tPm = new Point3d();
//                tPm.Interpolate(line.From, line.To, t[i]);
//                tempP.Add(tPm);
//            }
//            return tempP;
//        }

//        private static double[] LinSpace(double x1, double x2, int n)
//        {
//            //Generate a 1-D array of linearly spaced values
//            double step = (x2 - x1) / (n - 1);
//            double[] y = new double[n];
//            for (int i = 0; i < n; i++)
//            {
//                y[i] = x1 + step * i;
//            }
//            return y;
//        }

//        private void Shapefunctions(double L, double x, out Matrix<double> N, out Matrix<double> dN)
//        {
//            double N1 = -1 / L * (x - L);
//            double N2 = x / L;
//            double N3 = 1 - 3 * Math.Pow(x, 2) / Math.Pow(L, 2) + 2 * Math.Pow(x, 3) / Math.Pow(L, 3);
//            double N4 = x * (1 - 2 * x / L + Math.Pow(x, 2) / Math.Pow(L, 2));
//            double N5 = Math.Pow(x, 2) / Math.Pow(L, 2) * (3 - 2 * x / L);
//            double N6 = Math.Pow(x, 2) / L * (x / L - 1);

//            N = Matrix<double>.Build.DenseOfArray(new double[,] {
//                { N1, 0, 0,  0,  0,  0, N2, 0,  0,  0,  0,  0},
//                { 0, N3, 0,  0,  0, N4, 0, N5, 0,  0,  0, N6 },
//                { 0, 0, N3, 0, -N4, 0, 0, 0, N5, 0, -N6, 0},
//                { 0, 0, 0, N1, 0, 0, 0, 0, 0, N2, 0, 0} });

//            double dN1 = -1 / L;
//            double dN2 = 1 / L;
//            double dN3 = -6 * x / Math.Pow(L, 2) + 6 * Math.Pow(x, 2) / Math.Pow(L, 3);
//            double dN4 = 1 - 4 * x / L + 3 * Math.Pow(x, 2) / Math.Pow(L, 2);
//            double dN5 = 6 * x / Math.Pow(L, 2) - 6 * Math.Pow(x, 2) / Math.Pow(L, 3);
//            double dN6 = 3 * Math.Pow(x, 2) / Math.Pow(L, 2) - 2 * x / L;

//            dN = Matrix<double>.Build.DenseOfArray(new double[,] {
//            { dN1, 0, 0, 0, 0, 0, dN2, 0, 0, 0, 0, 0},
//            { 0, dN3, 0, 0, 0, dN4, 0, dN5, 0, 0, 0, dN6 },
//            { 0, 0, dN3, 0, -dN4, 0, 0, 0, dN5, 0, -dN6, 0},
//            { 0, 0, 0, dN1, 0, 0, 0, 0, 0, dN2, 0, 0} });
//        }

//        private void CalculateInternalStrainsAndStresses(Vector<double> def, List<Point3d> points, double E, List<Line> geometry, out List<double> internalStresses, out List<double> internalStrains)
//        {
//            //preallocating lists
//            internalStresses = new List<double>(geometry.Count);
//            internalStrains = new List<double>(geometry.Count);

//            foreach (Line line in geometry)
//            {
//                int index1 = points.IndexOf(new Point3d(Math.Round(line.From.X, 2), Math.Round(line.From.Y, 2), Math.Round(line.From.Z, 2)));
//                int index2 = points.IndexOf(new Point3d(Math.Round(line.To.X, 2), Math.Round(line.To.Y, 2), Math.Round(line.To.Z, 2)));

//                //fetching deformation of point
//                double x1 = def[index1 * 3 + 0];
//                double y1 = def[index1 * 3 + 1];
//                double z1 = def[index1 * 3 + 2];
//                double x2 = def[index2 * 3 + 0];
//                double y2 = def[index2 * 3 + 1];
//                double z2 = def[index2 * 3 + 2];

//                //new node coordinates for deformed nodes
//                double nx1 = points[index1].X + x1;
//                double ny1 = points[index1].X + y1;
//                double nz1 = points[index1].Z + z1;
//                double nx2 = points[index2].X + x2;
//                double ny2 = points[index2].X + y2;
//                double nz2 = points[index2].Z + z2;

//                //calculating dL = length of deformed line - original length of line
//                double dL = Math.Sqrt(Math.Pow((nx2 - nx1), 2) + Math.Pow((ny2 - ny1), 2) + Math.Pow((nz2 - nz1), 2)) - line.Length;

//                //calculating strain and stress
//                internalStrains.Add(dL / line.Length);
//                internalStresses.Add(internalStrains[internalStrains.Count - 1] * E);
//            }
//        }

//        private Vector<double> RestoreTotalDeformationVector(Vector<double> deformations_red, List<int> bdc_value)
//        {
//            Vector<double> def = Vector<double>.Build.Dense(bdc_value.Count);
//            for (int i = 0, j = 0; i < bdc_value.Count; i++)
//            {
//                if (bdc_value[i] == 1)
//                {
//                    def[i] = deformations_red[j];
//                    j++;
//                }
//            }
//            return def;
//        }

//        private void CreateReducedGlobalStiffnessMatrix(List<int> bdc_value, Matrix<double> K, List<double> load, out Matrix<double> K_red, out Vector<double> load_red)
//        {
//            K_red = Matrix<double>.Build.DenseOfMatrix(K);
//            List<double> load_redu = new List<double>(load);
//            for (int i = 0, j = 0; i < load.Count; i++)
//            {
//                if (bdc_value[i] == 0)
//                {
//                    K_red = K_red.RemoveRow(i - j);
//                    K_red = K_red.RemoveColumn(i - j);
//                    load_redu.RemoveAt(i - j);
//                    j++;
//                }
//            }
//            load_red = Vector<double>.Build.DenseOfEnumerable(load_redu);
//        }

//        private static bool IsDiagonalPositive(Matrix<double> A)
//        {
//            bool isPositive = true;
//            for (int i = 0; i < A.RowCount; i++)
//            {
//                if (A[i, i] <= 0)
//                {
//                    return isPositive = false;
//                }
//            }
//            return isPositive;
//        }

//        private void CheckSolvers(Matrix<double> K_red, Vector<double> load_red, int decimals, out string output_time)
//        {
//            string time;
//            //Unrounded
//            output_time = "Unrounded K_red" + Environment.NewLine;
//            bool issym = K_red.IsSymmetric();
//            bool diagposi = IsDiagonalPositive(K_red);
//            bool ishermitian = K_red.IsHermitian();
//            output_time += "IsSymmetric? " + issym + ", IsDiagonalPositive? " + diagposi + ", IsHermitian? " + ishermitian + Environment.NewLine + Environment.NewLine;

//            TrySolve(K_red, load_red, out time);
//            output_time += "Dense, unrounded K_red: " + Environment.NewLine + time + Environment.NewLine;

//            K_red = Matrix<double>.Build.SparseOfMatrix(K_red);
//            TrySolve(K_red, load_red, out time);
//            output_time += "Sparse, unrounded K_red: " + Environment.NewLine + time + Environment.NewLine;

//            //Rounded
//            K_red = Matrix<double>.Build.DenseOfMatrix(K_red);
//            K_red = RoundMatrix(K_red, decimals);
//            output_time += "Rounded to " + decimals + " decimals." + Environment.NewLine;

//            issym = K_red.IsSymmetric();
//            diagposi = IsDiagonalPositive(K_red);
//            ishermitian = K_red.IsHermitian();
//            output_time += "IsSymmetric? " + issym + ", IsDiagonalPositive? " + diagposi + ", IsHermitian? " + ishermitian + Environment.NewLine + Environment.NewLine;

//            TrySolve(K_red, load_red, out time);
//            output_time += "Dense, rounded K_red: " + Environment.NewLine + time + Environment.NewLine;

//            K_red = Matrix<double>.Build.SparseOfMatrix(K_red);
//            TrySolve(K_red, load_red, out time);
//            output_time += "Sparse, rounded K_red: " + Environment.NewLine + time + Environment.NewLine;
//            output_time += "=================END OF TEST=================" + Environment.NewLine + Environment.NewLine;
//        }

//        private void TrySolve(Matrix<double> A, Vector<double> load_red, out string time)
//        {
//            time = "TrySolve Start:" + Environment.NewLine;
//            long timer = 0;
//            Stopwatch watch = new Stopwatch();
//            try
//            {
//                watch.Start();
//                A.Solve(load_red);
//                watch.Stop();
//                timer = watch.ElapsedMilliseconds;
//                time += "Regular solve: " + timer.ToString() + Environment.NewLine;
//            }
//            catch (Exception)
//            {
//                time += "Regular solve raised exception" + Environment.NewLine;
//            }

//            if (A.GetType().Name == "SparseMatrix")
//            {
//                time += "End of TrySolve since Matrix is Sparse (other solvers are unsupported)" + Environment.NewLine;
//                return;
//            }

//            try
//            {
//                watch.Start();
//                A.Cholesky().Solve(load_red);
//                watch.Stop();
//                timer = watch.ElapsedMilliseconds - timer;
//                time += "Cholesky solve: " + timer.ToString() + Environment.NewLine;
//            }
//            catch (Exception)
//            {
//                time += "Cholesky solve raised exception" + Environment.NewLine;
//            }
//            try
//            {
//                watch.Start();
//                A.QR().Solve(load_red);
//                watch.Stop();
//                timer = watch.ElapsedMilliseconds - timer;
//                time += "QR solve: " + timer.ToString() + Environment.NewLine;
//            }
//            catch (Exception)
//            {
//                time += "QR solve raised exception" + Environment.NewLine;
//            }

//            try
//            {
//                watch.Start();
//                A.Svd().Solve(load_red);
//                watch.Stop();
//                timer = watch.ElapsedMilliseconds - timer;
//                time += "Svd solve: " + timer.ToString() + Environment.NewLine;
//            }
//            catch (Exception)
//            {
//                time += "Svd solve raised exception" + Environment.NewLine;
//            }
//            try
//            {
//                watch.Start();
//                A.LU().Solve(load_red);
//                watch.Stop();
//                timer = watch.ElapsedMilliseconds - timer;
//                time += "LU solve: " + timer.ToString() + Environment.NewLine;
//            }
//            catch (Exception)
//            {
//                time += "LU solve raised exception" + Environment.NewLine;
//            }
//        }

//        private string PrintMatrix(object B, string head)
//        {
//            Debug.WriteLine(head);
//            string ss = "";
//            if (B.GetType() == typeof(SparseMatrix))
//            {
//                Matrix<double> A = (Matrix<double>)B;
//                ss += A.RowCount.ToString() + "x" + A.ColumnCount.ToString() + "-matrix" + Environment.NewLine;
//                for (int i = 0; i < A.RowCount; i++)
//                {
//                    for (int j = 0; j < A.ColumnCount; j++)
//                    {
//                        ss += String.Format(" {0,15:0.0} ", A[i, j]);
//                    }
//                    ss += Environment.NewLine + Environment.NewLine + Environment.NewLine;
//                }
//                Debug.WriteLine(ss);
//                return ss;
//            }
//            else if (B.GetType() == typeof(DenseMatrix))
//            {
//                Matrix<double> A = (Matrix<double>)B;
//                ss += A.RowCount.ToString() + "x" + A.ColumnCount.ToString() + "-matrix" + Environment.NewLine;
//                for (int i = 0; i < A.RowCount; i++)
//                {
//                    for (int j = 0; j < A.ColumnCount; j++)
//                    {
//                        ss += String.Format(" {0,15:0.0} ", A[i, j]);
//                    }
//                    ss += Environment.NewLine + Environment.NewLine + Environment.NewLine;
//                }
//                Debug.WriteLine(ss);
//                return ss;
//            }
//            else if (B.GetType() == typeof(MathNet.Numerics.LinearAlgebra.Double.DenseVector))
//            {
//                Vector<double> A = (Vector<double>)B;
//                ss += "1x" + A.Count.ToString() + "-vector" + Environment.NewLine;
//                for (int i = 0; i < A.Count; i++)
//                {
//                    ss += String.Format(" {0,10:0.0} ", A[i]);
//                }
//                System.Diagnostics.Debug.WriteLine(ss);
//                return ss;
//            }
//            else if (B.GetType() == typeof(List<int>))
//            {
//                List<int> A = (List<int>)B;
//                ss += "1x" + A.Count.ToString() + "-list" + Environment.NewLine;
//                for (int i = 0; i < A.Count; i++)
//                {
//                    ss += String.Format(" {0,10:0.0} ", A[i]);
//                }
//                System.Diagnostics.Debug.WriteLine(ss);
//                return ss;
//            }
//            else if (B.GetType() == typeof(List<double>))
//            {
//                List<double> A = (List<double>)B;
//                ss += "1x" + A.Count.ToString() + "-list" + Environment.NewLine;
//                for (int i = 0; i < A.Count; i++)
//                {
//                    ss += String.Format(" {0,10:0.0} ", A[i]);
//                }
//                System.Diagnostics.Debug.WriteLine(ss);
//                return ss;
//            }
//            return ss;
//        }

//        private static Matrix<double> RoundMatrix(Matrix<double> A, int decs)
//        {
//            for (int i = 0; i < A.RowCount; i++)
//            {
//                for (int j = 0; j < A.ColumnCount; j++)
//                {
//                    A[i, j] = Math.Round(A[i, j], decs);
//                }
//            }
//            return A;
//        }

//        private void ElementStiffnessMatrix(Line currentLine, double E, double A, double Iy, double Iz, double J, double G, out Point3d p1, out Point3d p2, out Matrix<double> K_elem)
//        {
//            double L = Math.Round(currentLine.Length, 6);

//            p1 = new Point3d(Math.Round(currentLine.From.X, 2), Math.Round(currentLine.From.Y, 2), Math.Round(currentLine.From.Z, 2));
//            p2 = new Point3d(Math.Round(currentLine.To.X, 2), Math.Round(currentLine.To.Y, 2), Math.Round(currentLine.To.Z, 2));

//            double alpha = 0;

//            double cx = (p2.X - p1.X) / L;
//            double cy = (p2.Y - p1.Y) / L;
//            double cz = (p2.Z - p1.Z) / L;
//            double c1 = Math.Cos(alpha);
//            double s1 = Math.Sin(alpha);
//            double cxz = Math.Round(Math.Sqrt(Math.Pow(cx, 2) + Math.Pow(cz, 2)), 6);

//            Matrix<double> t;

//            if (Math.Round(cx, 6) == 0 && Math.Round(cz, 6) == 0)
//            {
//                t = Matrix<double>.Build.DenseOfArray(new double[,]
//            {
//                    {      0, cy,  0},
//                    { -cy*c1,  0, s1},
//                    {  cy*s1,  0, c1},
//            });
//            }
//            else
//            {
//                t = Matrix<double>.Build.DenseOfArray(new double[,]
//            {
//                    {                     cx,       cy,                   cz},
//                    {(-cx*cy*c1 - cz*s1)/cxz,   cxz*c1,(-cy*cz*c1+cx*s1)/cxz},
//                    {   (cx*cy*s1-cz*c1)/cxz,  -cxz*s1, (cy*cz*s1+cx*c1)/cxz},
//            });
//            }

//            var T = t.DiagonalStack(t);
//            T = T.DiagonalStack(T);

//            Matrix<double> T_T = T.Transpose();

//            double A1 = (E * A) / (L);

//            double kz1 = (12 * E * Iz) / (L * L * L);
//            double kz2 = (6 * E * Iz) / (L * L);
//            double kz3 = (4 * E * Iz) / L;
//            double kz4 = (2 * E * Iz) / L;

//            double ky1 = (12 * E * Iy) / (L * L * L);
//            double ky2 = (6 * E * Iy) / (L * L);
//            double ky3 = (4 * E * Iy) / L;
//            double ky4 = (2 * E * Iy) / L;

//            double C1 = (G * J) / L;

//            K_elem = DenseMatrix.OfArray(new double[,]
//            {
//                    { A1,    0,    0,    0,    0,    0,  -A1,    0,    0,    0,    0,    0 },
//                    {  0,  kz1,    0,    0,    0,  kz2,    0, -kz1,    0,    0,    0,  kz2 },
//                    {  0,    0,  ky1,    0, -ky2,    0,    0,    0, -ky1,    0, -ky2,    0 },
//                    {  0,    0,    0,   C1,    0,    0,    0,    0,    0,  -C1,    0,    0 },
//                    {  0,    0, -ky2,    0,  ky3,    0,    0,    0,  ky2,    0,  ky4,    0 },
//                    {  0,  kz2,    0,    0,    0,  kz3,    0, -kz2,    0,    0,    0,  kz4 },
//                    {-A1,    0,    0,    0,    0,    0,   A1,    0,    0,    0,    0,    0 },
//                    {  0, -kz1,    0,    0,    0, -kz2,    0,  kz1,    0,    0,    0, -kz2 },
//                    {  0,    0, -ky1,    0,  ky2,    0,    0,    0,  ky1,    0,  ky2,    0 },
//                    {  0,    0,    0,  -C1,    0,    0,    0,    0,    0,   C1,    0,    0 },
//                    {  0,    0, -ky2,    0,  ky4,    0,    0,    0,  ky2,    0,  ky3,    0 },
//                    {  0,  kz2,    0,    0,    0,  kz4,    0, -kz2,    0,    0,    0,  kz3 },
//            });

//            K_elem = K_elem.Multiply(T);
//            K_elem = T_T.Multiply(K_elem);
//        }

//        private Matrix<double> GlobalStiffnessMatrix(List<Line> geometry, List<Point3d> points, double E, double A, double Iy, double Iz, double J, double G)
//        {
//            int gdofs = points.Count * 6;
//            Matrix<double> K_tot = DenseMatrix.OfArray(new double[gdofs, gdofs]);

//            foreach (Line currentLine in geometry)
//            {
//                Matrix<double> K_elem;
//                Point3d p1;
//                Point3d p2;
//                ElementStiffnessMatrix(currentLine, E, A, Iy, Iz, J, G, out p1, out p2, out K_elem);

//                int node1 = points.IndexOf(p1);
//                int node2 = points.IndexOf(p2);

//                //Inputting values to correct entries in Global Stiffness Matrix
//                for (int i = 0; i < K_elem.RowCount / 2; i++)
//                {
//                    for (int j = 0; j < K_elem.ColumnCount / 2; j++)
//                    {
//                        //top left 3x3 of k-element matrix
//                        K_tot[node1 * 6 + i, node1 * 6 + j] += K_elem[i, j];
//                        //top right 3x3 of k-element matrix  
//                        K_tot[node1 * 6 + i, node2 * 6 + j] += K_elem[i, j + 6];
//                        //bottom left 3x3 of k-element matrix
//                        K_tot[node2 * 6 + i, node1 * 6 + j] += K_elem[i + 6, j];
//                        //bottom right 3x3 of k-element matrix
//                        K_tot[node2 * 6 + i, node2 * 6 + j] += K_elem[i + 6, j + 6];
//                    }
//                }
//            }
//            return K_tot;
//        }

//        private List<double> CreateLoadList(List<string> loadtxt, List<string> momenttxt, List<Point3d> points)
//        {
//            List<double> loads = new List<double>(new double[points.Count * 6]);
//            List<double> inputLoads = new List<double>();
//            List<Point3d> coordlist = new List<Point3d>();

//            for (int i = 0; i < loadtxt.Count; i++)
//            {
//                string coordstr = (loadtxt[i].Split(':')[0]);
//                string loadstr = (loadtxt[i].Split(':')[1]);

//                string[] coordstr1 = (coordstr.Split(','));
//                string[] loadstr1 = (loadstr.Split(','));

//                inputLoads.Add(Math.Round(double.Parse(loadstr1[0]), 2));
//                inputLoads.Add(Math.Round(double.Parse(loadstr1[1]), 2));
//                inputLoads.Add(Math.Round(double.Parse(loadstr1[2]), 2));

//                coordlist.Add(new Point3d(Math.Round(double.Parse(coordstr1[0]), 2), Math.Round(double.Parse(coordstr1[1]), 2), Math.Round(double.Parse(coordstr1[2]), 2)));
//            }

//            foreach (Point3d point in coordlist)
//            {
//                int i = points.IndexOf(point);
//                int j = coordlist.IndexOf(point);
//                loads[i * 6 + 0] = inputLoads[j * 3 + 0]; //is loads out of range? (doesn't seem to have been initialized with size yet)
//                loads[i * 6 + 1] = inputLoads[j * 3 + 1];
//                loads[i * 6 + 2] = inputLoads[j * 3 + 2];
//            }
//            inputLoads.Clear();
//            coordlist.Clear();
//            for (int i = 0; i < momenttxt.Count; i++) if (momenttxt[0] != "")
//                {
//                    string coordstr = (momenttxt[i].Split(':')[0]);
//                    string loadstr = (momenttxt[i].Split(':')[1]);

//                    string[] coordstr1 = (coordstr.Split(','));
//                    string[] loadstr1 = (loadstr.Split(','));

//                    inputLoads.Add(Math.Round(double.Parse(loadstr1[0]), 2));
//                    inputLoads.Add(Math.Round(double.Parse(loadstr1[1]), 2));
//                    inputLoads.Add(Math.Round(double.Parse(loadstr1[2]), 2));


//                    coordlist.Add(new Point3d(Math.Round(double.Parse(coordstr1[0]), 2), Math.Round(double.Parse(coordstr1[1]), 2), Math.Round(double.Parse(coordstr1[2]), 2)));
//                }

//            foreach (Point3d point in coordlist)
//            {
//                int i = points.IndexOf(point);
//                int j = coordlist.IndexOf(point);
//                loads[i * 6 + 3] = inputLoads[j * 3 + 0];
//                loads[i * 6 + 4] = inputLoads[j * 3 + 1];
//                loads[i * 6 + 5] = inputLoads[j * 3 + 2];
//            }
//            return loads;
//        }

//        private List<int> CreateBDCList(List<string> bdctxt, List<Point3d> points)
//        {
//            List<int> bdc_value = new List<int>(new int[points.Count * 6]);
//            List<int> bdcs = new List<int>();
//            List<Point3d> bdc_points = new List<Point3d>(); //Coordinates relating til bdc_value in for (eg. x y z)

//            //Parse string input
//            for (int i = 0; i < bdctxt.Count; i++)
//            {
//                string coordstr = (bdctxt[i].Split(':')[0]);
//                string bdcstr = (bdctxt[i].Split(':')[1]);

//                string[] coordstr1 = (coordstr.Split(','));
//                string[] bdcstr1 = (bdcstr.Split(','));

//                bdc_points.Add(new Point3d(Math.Round(double.Parse(coordstr1[0]), 2), Math.Round(double.Parse(coordstr1[1]), 2), Math.Round(double.Parse(coordstr1[2]), 2)));

//                bdcs.Add(int.Parse(bdcstr1[0]));
//                bdcs.Add(int.Parse(bdcstr1[1]));
//                bdcs.Add(int.Parse(bdcstr1[2]));
//                bdcs.Add(int.Parse(bdcstr1[3]));
//                bdcs.Add(int.Parse(bdcstr1[4]));
//                bdcs.Add(int.Parse(bdcstr1[5]));
//            }


//            //Format to correct entries in bdc_value
//            for (int i = 0; i < points.Count; i++)
//            {
//                Point3d tempP = points[i];

//                if (bdc_points.Contains(tempP))
//                {
//                    bdc_value[i * 6 + 0] = bdcs[bdc_points.IndexOf(tempP) * 6 + 0];
//                    bdc_value[i * 6 + 1] = bdcs[bdc_points.IndexOf(tempP) * 6 + 1];
//                    bdc_value[i * 6 + 2] = bdcs[bdc_points.IndexOf(tempP) * 6 + 2];
//                    bdc_value[i * 6 + 3] = bdcs[bdc_points.IndexOf(tempP) * 6 + 3];
//                    bdc_value[i * 6 + 4] = bdcs[bdc_points.IndexOf(tempP) * 6 + 4];
//                    bdc_value[i * 6 + 5] = bdcs[bdc_points.IndexOf(tempP) * 6 + 5];
//                }
//                else
//                {
//                    bdc_value[i * 6 + 0] = 1;
//                    bdc_value[i * 6 + 1] = 1;
//                    bdc_value[i * 6 + 2] = 1;
//                    bdc_value[i * 6 + 3] = 1;
//                    bdc_value[i * 6 + 4] = 1;
//                    bdc_value[i * 6 + 5] = 1;
//                }
//            }
//            return bdc_value;
//        }

//        private void SetMaterial(string mattxt, out double E, out double A, out double Iy, out double Iz, out double J, out double G)
//        {
//            string[] matProp = (mattxt.Split(','));

//            E = (Math.Round(double.Parse(matProp[0]), 2));
//            A = (Math.Round(double.Parse(matProp[1]), 2));
//            Iy = (Math.Round(double.Parse(matProp[2]), 2));
//            Iz = (Math.Round(double.Parse(matProp[3]), 2));
//            G = (Math.Round(double.Parse(matProp[4]), 2));
//            J = Iy + Iz;
//        }

//        private List<Point3d> CreatePointList(List<Line> geometry)
//        {
//            List<Point3d> points = new List<Point3d>();
//            foreach (Line line in geometry) //adds point unless it already exists in pointlist
//            {
//                Point3d tempFrom = new Point3d(Math.Round(line.From.X, 2), Math.Round(line.From.Y, 2), Math.Round(line.From.Z, 2));
//                Point3d tempTo = new Point3d(Math.Round(line.To.X, 2), Math.Round(line.To.Y, 2), Math.Round(line.To.Z, 2));

//                if (!points.Contains(tempFrom))
//                {
//                    points.Add(tempFrom);
//                }
//                if (!points.Contains(tempTo))
//                {
//                    points.Add(tempTo);
//                }
//            }
//            return points;
//        }

//        protected override System.Drawing.Bitmap Icon
//        {
//            get
//            {
//                return Properties.Resources.Calc;
//            }
//        }

//        public override Guid ComponentGuid
//        {
//            get { return new Guid("d636ebc9-0d19-44d5-a3ad-cec704b82323"); }
//        }


//        /// Component Visual//
//        public class Attributes_Custom : Grasshopper.Kernel.Attributes.GH_ComponentAttributes
//        {
//            public Attributes_Custom(GH_Component owner) : base(owner) { }
//            protected override void Layout()
//            {
//                base.Layout();

//                Rectangle rec0 = GH_Convert.ToRectangle(Bounds);

//                rec0.Height += 22;

//                Rectangle rec1 = rec0;
//                rec1.X = rec0.Left + 1;
//                rec1.Y = rec0.Bottom - 22;
//                rec1.Width = (rec0.Width) / 3 + 1;
//                rec1.Height = 22;
//                rec1.Inflate(-2, -2);

//                Rectangle rec2 = rec1;
//                rec2.X = rec1.Right + 2;

//                Bounds = rec0;
//                ButtonBounds = rec1;
//                ButtonBounds2 = rec2;

//            }

//            GH_Palette xColor = GH_Palette.Black;
//            GH_Palette yColor = GH_Palette.Grey;

//            private Rectangle ButtonBounds { get; set; }
//            private Rectangle ButtonBounds2 { get; set; }
//            private Rectangle ButtonBounds3 { get; set; }

//            protected override void Render(GH_Canvas canvas, Graphics graphics, GH_CanvasChannel channel)
//            {
//                base.Render(canvas, graphics, channel);
//                if (channel == GH_CanvasChannel.Objects)
//                {
//                    GH_Capsule button = GH_Capsule.CreateTextCapsule(ButtonBounds, ButtonBounds, xColor, "Run", 3, 0);
//                    button.Render(graphics, Selected, false, false);
//                    button.Dispose();
//                }
//                if (channel == GH_CanvasChannel.Objects)
//                {
//                    GH_Capsule button2 = GH_Capsule.CreateTextCapsule(ButtonBounds2, ButtonBounds2, yColor, "Run Test", 2, 0);
//                    button2.Render(graphics, Selected, Owner.Locked, false);
//                    button2.Dispose();
//                }
//            }

//            public override GH_ObjectResponse RespondToMouseDown(GH_Canvas sender, GH_CanvasMouseEvent e)
//            {
//                if (e.Button == MouseButtons.Left)
//                {
//                    RectangleF rec = ButtonBounds;
//                    if (rec.Contains(e.CanvasLocation))
//                    {
//                        switchColor("Run");
//                        if (xColor == GH_Palette.Black) { CalcComponent.setStart("Run", true); }
//                        if (xColor == GH_Palette.Grey) { CalcComponent.setStart("Run", false); }
//                        sender.Refresh();
//                        return GH_ObjectResponse.Handled;
//                    }
//                    rec = ButtonBounds2;
//                    if (rec.Contains(e.CanvasLocation))
//                    {
//                        switchColor("Run Test");
//                        if (yColor == GH_Palette.Black) { CalcComponent.setStart("Run Test", true); }
//                        if (yColor == GH_Palette.Grey) { CalcComponent.setStart("Run Test", false); }
//                        sender.Refresh();
//                        return GH_ObjectResponse.Handled;
//                    }
//                }
//                return base.RespondToMouseDown(sender, e);
//            }

//            private void switchColor(string button)
//            {
//                if (button == "Run")
//                {
//                    if (xColor == GH_Palette.Black) { xColor = GH_Palette.Grey; }
//                    else { xColor = GH_Palette.Black; }
//                }
//                else if (button == "Run Test")
//                {
//                    if (yColor == GH_Palette.Black) { yColor = GH_Palette.Grey; }
//                    else { yColor = GH_Palette.Black; }
//                }
//            }
//        }
//    }
//}
#endregion