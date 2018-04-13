using System;
using System.Collections.Generic;

using Grasshopper.Kernel;
using Rhino.Geometry;

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

        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddLineParameter("Lines", "LNS", "Geometry, in form of Lines)", GH_ParamAccess.list);
            pManager.AddTextParameter("Boundary Conditions", "BDC", "Boundary Conditions in form x,y,z,vx,vy,vz,rx,ry,rz", GH_ParamAccess.list);
            pManager.AddTextParameter("Material properties", "Mat", "Material Properties", GH_ParamAccess.item, "210000,3600,4920000,4920000,79300");
            pManager.AddTextParameter("PointLoads", "PL", "Load given as Vector [N]", GH_ParamAccess.list);
            pManager.AddTextParameter("PointMoment", "PM", "Moment set in a point in [Nm]", GH_ParamAccess.list, "");
            pManager.AddBooleanParameter("Start calculations", "SC", "Set true to start calculations", GH_ParamAccess.item, false);
            pManager.AddBooleanParameter("Solver test?", "ST", "Set true to start solver performance test", GH_ParamAccess.item, false);
        }

        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddNumberParameter("Deformations", "Def", "Deformations", GH_ParamAccess.list);
            pManager.AddNumberParameter("Reactions", "R", "Reaction Forces", GH_ParamAccess.list);
            pManager.AddNumberParameter("Element stresses", "Strs", "The Stress in each element", GH_ParamAccess.list);
            pManager.AddNumberParameter("Element strains", "Strn", "The Strain in each element", GH_ParamAccess.list);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            #region Fetch input
            //Expected inputs
            List<Line> geometry = new List<Line>();         //Initial Geometry of lines
            List<string> bdctxt = new List<string>();       //Boundary conditions in string format
            List<string> loadtxt = new List<string>();      //loads in string format
            List<string> momenttxt = new List<string>();    //Moments in string format
            string mattxt = "";                             //Material in string format
            bool startCalc = false;
            bool startTest = false;


            //Set expected inputs from Indata
            if (!DA.GetDataList(0, geometry)) return;       //sets geometry
            if (!DA.GetDataList(1, bdctxt)) return;         //sets boundary conditions as string
            if (!DA.GetData(2, ref mattxt)) return;         //sets material properties as string
            if (!DA.GetDataList(3, loadtxt)) return;        //sets load as string
            if (!DA.GetDataList(4, momenttxt)) return;      //sets moment as string
            if (!DA.GetData(5, ref startCalc)) return;
            if (!DA.GetData(6, ref startTest)) return;
            #endregion

            if (startCalc)
            {
                //Interpret and set material parameters
                double E;       //Material Young's modulus, initial value 210000 [MPa]
                double A;       //Area for each element in same order as geometry, initial value CFS100x100 3600 [mm^2]
                double Iy;      //Moment of inertia about local y axis, initial value 4.92E6 [mm^4]
                double Iz;      //Moment of inertia about local z axis, initial value 4.92E6 [mm^4]
                double J;       //Polar moment of inertia
                double G;       //Shear modulus, initial value 79300 [mm^4]
                SetMaterial(mattxt, out E, out A, out Iy, out Iz, out J, out G);
                
                #region Prepares geometry, boundary conditions and loads for calculation
                //List all nodes (every node only once), numbering them according to list index
                List<Point3d> points = CreatePointList(geometry);


                //Interpret the BDC inputs (text) and create list of boundary condition (1/0 = free/clamped) for each dof.
                List<int> bdc_value = CreateBDCList(bdctxt, points);


                //Interpreting input load (text) and creating load list (do uble)
                List<double> load = CreateLoadList(loadtxt, momenttxt, points);
                #endregion

                #region Create global and reduced stiffness matrix
                //Create global stiffness matrix
                Matrix<double> K_tot = CreateGlobalStiffnessMatrix(geometry, points, E, A, Iy, Iz, J, G);

                //Create reduced K-matrix and reduced load list (removed free dofs)
                Matrix<double> K_red;
                Vector<double> load_red;
                CreateReducedGlobalStiffnessMatrix(bdc_value, K_tot, load, out K_red, out load_red);
                #endregion

                #region Solver Performance Test
                if (startTest)
                {
                    string output_time = "";
                    string performanceResult = "=================START OF TEST=================" + Environment.NewLine;
                    performanceResult += "Number of lines: " + geometry.Count.ToString() + Environment.NewLine;
                    string tester = "";

                    //checking for error with writing to file (skip test if unable to write)
                    try
                    {
                        using (System.IO.StreamWriter file =
                        new System.IO.StreamWriter(@"solverTest.txt", true))
                        {
                            file.WriteLine(tester);
                        }
                    }
                    //create new file if solverTest.txt does not exist
                    catch (System.IO.DirectoryNotFoundException)
                    {
                        System.IO.File.WriteAllText(@"solverTest.txt", tester);
                    }
                    //other exception (no write access?)
                    catch (Exception)
                    {
                        AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Write to file error! Create file solverTest.txt in Koala folder and/or Rhinoceros 5.exe folder");
                        return;
                    }

                    int decimals = 2;
                    CheckSolvers(K_red, load_red, decimals, out output_time);

                    performanceResult += output_time;

                    //append result to txt-file (at buildpath)
                    try
                    {
                        using (System.IO.StreamWriter file =
                        new System.IO.StreamWriter(@"solverTest.txt", true))
                        {
                            file.WriteLine(performanceResult);
                        }
                    }
                    //create new file if solverTest.txt does not exist
                    catch (System.IO.DirectoryNotFoundException)
                    {
                        System.IO.File.WriteAllText(@"\solverTest.txt", output_time);
                    }
                    
                }
                #endregion
                
                #region Calculate deformations, reaction forces and internal strains and stresses
                //Calculate deformations
                Vector<double> def_reduced;
                def_reduced = K_red.Cholesky().Solve(load_red);

                //Add the clamped dofs (= 0) to the deformations list
                Vector<double> def_tot = RestoreTotalDeformationVector(def_reduced, bdc_value);

                //Calculate the reaction forces from the deformations
                Vector<double> reactions = K_tot.Multiply(def_tot);
                
                //Calculate the internal strains and stresses in each member
                List<double> internalStresses;
                List<double> internalStrains;
                CalculateInternalStrainsAndStresses(def_tot, points, E, geometry, out internalStresses, out internalStrains);
                #endregion

                DA.SetDataList(0, def_tot);
                DA.SetDataList(1, reactions);
                DA.SetDataList(2, internalStresses);
                DA.SetDataList(3, internalStrains);
            }
        } //End of main component

        private void CalculateInternalStrainsAndStresses(Vector<double> def, List<Point3d> points, double E, List<Line> geometry, out List<double> internalStresses, out List<double> internalStrains)
        {
            //preallocating lists
            internalStresses = new List<double>(geometry.Count);
            internalStrains = new List<double>(geometry.Count);

            foreach (Line line in geometry)
            {
                int index1 = points.IndexOf(new Point3d(Math.Round(line.From.X, 2), Math.Round(line.From.Y, 2), Math.Round(line.From.Z, 2)));
                int index2 = points.IndexOf(new Point3d(Math.Round(line.To.X, 2), Math.Round(line.To.Y, 2), Math.Round(line.To.Z, 2)));

                //fetching deformation of point
                double x1 = def[index1 * 3 + 0];
                double y1 = def[index1 * 3 + 1];
                double z1 = def[index1 * 3 + 2];
                double x2 = def[index2 * 3 + 0];
                double y2 = def[index2 * 3 + 1];
                double z2 = def[index2 * 3 + 2];

                //new node coordinates for deformed nodes
                double nx1 = points[index1].X + x1;
                double ny1 = points[index1].X + y1;
                double nz1 = points[index1].Z + z1;
                double nx2 = points[index2].X + x2;
                double ny2 = points[index2].X + y2;
                double nz2 = points[index2].Z + z2;

                //calculating dL = length of deformed line - original length of line
                double dL = Math.Sqrt(Math.Pow((nx2 - nx1), 2) + Math.Pow((ny2 - ny1), 2) + Math.Pow((nz2 - nz1), 2)) - line.Length;

                //calculating strain and stress
                internalStrains.Add(dL / line.Length);
                internalStresses.Add(internalStrains[internalStrains.Count - 1] * E);
            }
        }

        private Vector<double> RestoreTotalDeformationVector(Vector<double> deformations_red, List<int> bdc_value)
        {
            Vector<double> def = Vector<double>.Build.Dense(bdc_value.Count);
            for (int i = 0, j = 0; i < bdc_value.Count; i++)
            {
                if (bdc_value[i] == 1)
                {
                    def[i] = deformations_red[j];
                    j++;
                }
            }
            return def;
        }

        private void CreateReducedGlobalStiffnessMatrix(List<int> bdc_value, Matrix<double> K, List<double> load, out Matrix<double> K_red, out Vector<double> load_red)
        {
            K_red = Matrix<double>.Build.DenseOfMatrix(K);
            List<double> load_redu = new List<double>(load);
            for (int i = 0, j = 0; i < load.Count; i++)
            {
                if (bdc_value[i] == 0)
                {
                    K_red = K_red.RemoveRow(i - j);
                    K_red = K_red.RemoveColumn(i - j);
                    load_redu.RemoveAt(i - j);
                    j++;
                }
            }
            load_red = Vector<double>.Build.DenseOfEnumerable(load_redu);
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

        private Matrix<double> CreateGlobalStiffnessMatrix(List<Line> geometry, List<Point3d> points, double E, double A, double Iy, double Iz, double J, double G)
        {
            int gdofs = points.Count * 6;
            Matrix<double> K_tot = DenseMatrix.OfArray(new double[gdofs, gdofs]);

            foreach (Line currentLine in geometry)
            {
                double L = Math.Round(currentLine.Length, 6);

                Point3d p1 = new Point3d(Math.Round(currentLine.From.X, 2), Math.Round(currentLine.From.Y, 2), Math.Round(currentLine.From.Z, 2));
                Point3d p2 = new Point3d(Math.Round(currentLine.To.X, 2), Math.Round(currentLine.To.Y, 2), Math.Round(currentLine.To.Z, 2));

                double alpha = 0;

                double cx = (p2.X - p1.X) / L;
                double cy = (p2.Y - p1.Y) / L;
                double cz = (p2.Z - p1.Z) / L;
                double c1 = Math.Cos(alpha);
                double s1 = Math.Sin(alpha);
                double cxz = Math.Round(Math.Sqrt(Math.Pow(cx, 2) + Math.Pow(cz, 2)), 6);

                Matrix<double> gamma;

                if (Math.Round(cx, 6) == 0 && Math.Round(cz, 6) == 0)
                {
                    gamma = Matrix<double>.Build.DenseOfArray(new double[,]
                {
                    {      0, cy,  0},
                    { -cy*c1,  0, s1},
                    {  cy*s1,  0, c1},
                });
                }
                else
                {
                    gamma = Matrix<double>.Build.DenseOfArray(new double[,]
                {
                    {                     cx,       cy,                   cz},
                    {(-cx*cy*c1 - cz*s1)/cxz,   cxz*c1,(-cy*cz*c1+cx*s1)/cxz},
                    {   (cx*cy*s1-cz*c1)/cxz,  -cxz*s1, (cy*cz*s1+cx*c1)/cxz},
                });
                }

                var bd = Matrix<double>.Build;

                Matrix<double> T1;
                T1 = gamma.Append(bd.Dense(3, 9));
                Matrix<double> T2;
                T2 = bd.Dense(3, 3).Append(gamma);
                T2 = T2.Append(bd.Dense(3, 6));
                Matrix<double> T3;
                T3 = bd.Dense(3, 6).Append(gamma);
                T3 = T3.Append(bd.Dense(3, 3));
                Matrix<double> T4;
                T4 = bd.Dense(3, 9).Append(gamma);
                Matrix<double> T;
                T = T1.Stack(T2);
                T = T.Stack(T3);
                T = T.Stack(T4);

                //Matrix<double> T = SparseMatrix.OfArray(new double[,]
                //{
                //    { cx, cy, cz, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
                //    { cx, cy, cz, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
                //    { cx, cy, cz, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
                //    { 0, 0, 0, cx, cy, cz, 0, 0, 0, 0, 0, 0 },
                //    { 0, 0, 0, cx, cy, cz, 0, 0, 0, 0, 0, 0 },
                //    { 0, 0, 0, cx, cy, cz, 0, 0, 0, 0, 0, 0 },
                //    { 0, 0, 0, 0, 0, 0, cx, cy, cz, 0, 0, 0 },
                //    { 0, 0, 0, 0, 0, 0, cx, cy, cz, 0, 0, 0 },
                //    { 0, 0, 0, 0, 0, 0, cx, cy, cz, 0, 0, 0 },
                //    { 0, 0, 0, 0, 0, 0, 0, 0, 0, cx, cy, cz },
                //    { 0, 0, 0, 0, 0, 0, 0, 0, 0, cx, cy, cz },
                //    { 0, 0, 0, 0, 0, 0, 0, 0, 0, cx, cy, cz },
                //});

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

                Matrix<double> K_elem = DenseMatrix.OfArray(new double[,]
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

                K_elem = K_elem.Multiply(T);
                K_elem = T_T.Multiply(K_elem);

                int node1 = points.IndexOf(p1);
                int node2 = points.IndexOf(p2);

                //System.Diagnostics.Debug.WriteLine("Node1: " + node1.ToString() + ", Node2: " + node2.ToString());

                //PrintMatrix(K_elem,"K_elem");

                //Inputting values to correct entries in Global Stiffness Matrix
                for (int i = 0; i < K_elem.RowCount / 2; i++)
                {
                    //top left 3x3 of k-element matrix
                    for (int j = 0; j < K_elem.ColumnCount / 2; j++)
                    {
                        K_tot[node1 * 6 + i, node1 * 6 + j] += K_elem[i, j];
                    }
                    //top right 3x3 of k-element matrix  
                    for (int j = 0; j < K_elem.ColumnCount / 2; j++)
                    {
                        K_tot[node1 * 6 + i, node2 * 6 + j] += K_elem[i, j + 6];
                    }
                    //bottom left 3x3 of k-element matrix
                    for (int j = 0; j < K_elem.ColumnCount / 2; j++)
                    {
                        K_tot[node2 * 6 + i, node1 * 6 + j] += K_elem[i + 6, j];
                    }
                    //bottom right 3x3 of k-element matrix
                    for (int j = 0; j < K_elem.ColumnCount / 2; j++)
                    {
                        K_tot[node2 * 6 + i, node2 * 6 + j] += K_elem[i + 6, j + 6];
                    }
                }
            }
            return K_tot;
        }

        private List<double> CreateLoadList(List<string> loadtxt, List<string> momenttxt, List<Point3d> points)
        {
            List<double> loads = new List<double>(new double[points.Count * 6]);
            List<double> inputLoads = new List<double>();
            List<Point3d> coordlist = new List<Point3d>();

            for (int i = 0; i < loadtxt.Count; i++)
            {
                string coordstr = (loadtxt[i].Split(':')[0]);
                string loadstr = (loadtxt[i].Split(':')[1]);

                string[] coordstr1 = (coordstr.Split(','));
                string[] loadstr1 = (loadstr.Split(','));

                inputLoads.Add(Math.Round(double.Parse(loadstr1[0]), 2));
                inputLoads.Add(Math.Round(double.Parse(loadstr1[1]), 2));
                inputLoads.Add(Math.Round(double.Parse(loadstr1[2]), 2));

                coordlist.Add(new Point3d(Math.Round(double.Parse(coordstr1[0]), 2), Math.Round(double.Parse(coordstr1[1]), 2), Math.Round(double.Parse(coordstr1[2]), 2)));
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

                    inputLoads.Add(Math.Round(double.Parse(loadstr1[0]), 2));
                    inputLoads.Add(Math.Round(double.Parse(loadstr1[1]), 2));
                    inputLoads.Add(Math.Round(double.Parse(loadstr1[2]), 2));


                    coordlist.Add(new Point3d(Math.Round(double.Parse(coordstr1[0]), 2), Math.Round(double.Parse(coordstr1[1]), 2), Math.Round(double.Parse(coordstr1[2]), 2)));
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

        private List<int> CreateBDCList(List<string> bdctxt, List<Point3d> points)
        {
            List<int> bdc_value = new List<int>(new int[points.Count * 6]);
            List<int> bdcs = new List<int>();
            List<Point3d> bdc_points = new List<Point3d>(); //Coordinates relating til bdc_value in for (eg. x y z)

            //Parse string input
            for (int i = 0; i < bdctxt.Count; i++)
            {
                string coordstr = (bdctxt[i].Split(':')[0]);
                string bdcstr = (bdctxt[i].Split(':')[1]);

                string[] coordstr1 = (coordstr.Split(','));
                string[] bdcstr1 = (bdcstr.Split(','));

                bdc_points.Add(new Point3d(Math.Round(double.Parse(coordstr1[0]), 2), Math.Round(double.Parse(coordstr1[1]), 2), Math.Round(double.Parse(coordstr1[2]), 2)));

                bdcs.Add(int.Parse(bdcstr1[0]));
                bdcs.Add(int.Parse(bdcstr1[1]));
                bdcs.Add(int.Parse(bdcstr1[2]));
                bdcs.Add(int.Parse(bdcstr1[3]));
                bdcs.Add(int.Parse(bdcstr1[4]));
                bdcs.Add(int.Parse(bdcstr1[5]));
            }


            //Format to correct entries in bdc_value
            for (int i = 0; i < points.Count; i++)
            {
                Point3d tempP = points[i];

                if (bdc_points.Contains(tempP))
                {
                    bdc_value[i * 6 + 0] = bdcs[bdc_points.IndexOf(tempP) * 6 + 0];
                    bdc_value[i * 6 + 1] = bdcs[bdc_points.IndexOf(tempP) * 6 + 1];
                    bdc_value[i * 6 + 2] = bdcs[bdc_points.IndexOf(tempP) * 6 + 2];
                    bdc_value[i * 6 + 3] = bdcs[bdc_points.IndexOf(tempP) * 6 + 3];
                    bdc_value[i * 6 + 4] = bdcs[bdc_points.IndexOf(tempP) * 6 + 4];
                    bdc_value[i * 6 + 5] = bdcs[bdc_points.IndexOf(tempP) * 6 + 5];
                }
                else
                {
                    bdc_value[i * 6 + 0] = 1;
                    bdc_value[i * 6 + 1] = 1;
                    bdc_value[i * 6 + 2] = 1;
                    bdc_value[i * 6 + 3] = 1;
                    bdc_value[i * 6 + 4] = 1;
                    bdc_value[i * 6 + 5] = 1;
                }
            }
            return bdc_value;
        }

        private void SetMaterial(string mattxt, out double E, out double A, out double Iy, out double Iz, out double J, out double G)
        {
            string[] matProp = (mattxt.Split(','));

            E = (Math.Round(double.Parse(matProp[0]), 2));
            A = (Math.Round(double.Parse(matProp[1]), 2));
            Iy = (Math.Round(double.Parse(matProp[2]), 2));
            Iz = (Math.Round(double.Parse(matProp[3]), 2));
            G = (Math.Round(double.Parse(matProp[4]), 2));
            J = Iy + Iz;
        }

        private List<Point3d> CreatePointList(List<Line> geometry)
        {
            List<Point3d> points = new List<Point3d>();
            foreach (Line line in geometry) //adds point unless it already exists in pointlist
            {
                Point3d tempFrom = new Point3d(Math.Round(line.From.X, 2), Math.Round(line.From.Y, 2), Math.Round(line.From.Z, 2));
                Point3d tempTo = new Point3d(Math.Round(line.To.X, 2), Math.Round(line.To.Y, 2), Math.Round(line.To.Z, 2));

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
    }
}