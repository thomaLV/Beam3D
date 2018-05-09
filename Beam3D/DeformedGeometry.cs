using System;
using System.Collections.Generic;

using Grasshopper.Kernel;
using Rhino.Geometry;
using System.Drawing;
using Grasshopper.GUI.Canvas;
using System.Windows.Forms;
using Grasshopper.GUI;

using MathNet.Numerics.LinearAlgebra;

namespace Beam3D
{
    public class DeformedGeometry : GH_Component
    {
        public DeformedGeometry()
          : base("DeformedGeometry", "DefG",
              "Description",
              "Koala", "3D Beam")
        {
        }

        //Initialize startcondition
        static bool startDef = true;
        static bool p2 = true;
        static bool p3 = false;


        //Method to allow C# hanging of variables via GUI (see Component Visual)
        public static void setToggles(string s, bool i)
        {
            if (s == "Run")
            {
                startDef = i;
            }
            else if (s == "2nd")
            {
                p2 = i;
            }
            else if (s == "3rd")
            {
                p3 = i;
            }
            Grasshopper.Instances.ActiveCanvas.Document.ExpireSolution();
            Grasshopper.Instances.ActiveCanvas.Document.NewSolution(false);
        }

        public override void CreateAttributes()
        {
            m_attributes = new Attributes_Custom(this);
        }

        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddNumberParameter("Deformation", "Def", "Deformations from 3DBeamCalc", GH_ParamAccess.list);
            pManager.AddLineParameter("Geometry", "G", "Input Geometry (Line format)", GH_ParamAccess.list);
            pManager.AddNumberParameter("Scale", "S", "The Scale Factor for Deformation", GH_ParamAccess.item, 1000);
            pManager.AddIntegerParameter("Elements", "n", "No. of Elements", GH_ParamAccess.item, 4);
        }

        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddCurveParameter("Deformed Geometry", "Def.G.", "Deformed Geometry as List of Lines", GH_ParamAccess.list);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            if (startDef)
            {
                #region Fetch input
                //Expected inputs and outputs
                List<double> def = new List<double>();          //input doubles of displacements/rotations
                List<Line> geometry = new List<Line>();         //input list of original node coordinates
                double scale = 1000;                            //input deformation scale (default == 1000)                
                int n = 4;                                      //input no. of elements (default == 4)

                List<Curve> defGeometry = new List<Curve>();    //output deformed geometry
                List<Point3d> defPoints = new List<Point3d>();  //output deformed element nodes (remove?)

                var m = Matrix<double>.Build;
                var v = Vector<double>.Build;

                //Set expected inputs from Indata
                if (!DA.GetDataList(0, def)) return;
                if (!DA.GetDataList(1, geometry)) return;
                if (!DA.GetData(2, ref scale)) return;
                if (!DA.GetData(3, ref n)) return;
                #endregion

                //List all nodes (every node only once), numbering them according to list index
                List<Point3d> points = CreatePointList(geometry);


                #region Create geometry (linear)
                if (n == 1)
                {
                    int index = 0;
                    //loops through all points and scales x-, y- and z-dir

                    //u(x) = Na, N = shape func, a = nodal values (dof) 
                    foreach (Point3d point in points)
                    {
                        //fetch global x,y,z placement of point
                        double x = point.X;
                        double y = point.Y;
                        double z = point.Z;

                        //scales x and z according to input Scale
                        defPoints.Add(new Point3d(x + scale * def[index], y + scale * def[index + 1], z + scale * def[index + 2]));
                        index += 6;
                    }

                    //creates deformed geometry based on initial geometry placement
                    foreach (Line line in geometry)
                    {
                        //fetches index of original start and endpoint
                        int i1 = points.IndexOf(line.From);
                        int i2 = points.IndexOf(line.To);

                        //creates new line based on scaled deformation of said points
                        var tempL = new Line(defPoints[i1], defPoints[i2]);
                        defGeometry.Add(tempL.ToNurbsCurve());
                    }


                    //Set output data
                    DA.SetDataList(0, defGeometry);
                    return; //return to grasshopper since solved linearly
                }
                #endregion

                #region Create geometry
                Matrix<double> N, B;
                foreach (Line line in geometry)
                {
                    //fetches index of original start and endpoint
                    int i1 = points.IndexOf(line.From);
                    int i2 = points.IndexOf(line.To);

                    //create 12x1 deformation vector for element (6dofs), scaled and populated with existing deformations
                    var u = Vector<double>.Build.Dense(12);
                    for (int j = 0; j < 6; j++)
                    {
                        u[j] = def[i1*6 + j];
                        u[j + 6] = def[i2*6 + j];
                    }
                    u = scale * u;

                    //interpolate points between line.From and line.To
                    List<Point3d> tempP;
                    InterpolatePoints(line, n, out tempP);


                    double L = points[i1].DistanceTo(points[i2]);   //L is distance from startnode to endnode
                    var x = v.Dense(n + 1);                       //maybe this should be projected x instead???
                    for (int j = 0; j < n + 1; j++)
                    {
                        x[j] = j * L / n;
                    }


                    //Calculate 6 dofs for all new elements using shape functions (n+1 elements)
                    Matrix<double> disp = m.Dense(n + 1, 4);
                    Matrix<double> rot = m.Dense(n + 1, 4);

                    for (int j = 0; j < n + 1; j++)          //x are points inbetween (?)
                    {
                        Shapefunctions(L, x[j], out N, out B);
                        disp.SetRow(j, N.Multiply(u));
                        rot.SetRow(j, B.Multiply(u));
                    }

                    //Calculate new nodal points
                    for (int j = 0; j < n + 1; j++)
                    {
                        //original xyz                        
                        var tP = tempP[j];

                        //add displacement
                        //tP.X += disp[j, 0];
                        //tP.Y += disp[j, 1];
                        //tP.Z += disp[j, 2];

                        tP.X = tP.X + disp[j, 0] + tP.Z * Math.Cos(Math.PI / 2 - rot[j, 2]) + tP.Y * Math.Cos(Math.PI / 2 - rot[j, 1]);
                        tP.Y = -Math.Cos(disp[j, 3]) * tP.Y * Math.Sin(Math.PI / 2 - rot[j, 1]) + Math.Sin(disp[j, 3]) * tP.Z + disp[j, 1];
                        //tP.Z = disp[j, 2] + tP.Z * Math.Sin(rot[j, 2]);
                        tP.Z = -Math.Sin(disp[j, 3]) * tP.Y + Math.Cos(disp[j, 3]) * tP.Z * Math.Sin(Math.PI / 2 - rot[j, 3]) + disp[j, 2];

                        //replace previous xyz with displaced xyz
                        tempP[j] = tP;
                    }

                    //Create Nurbscurve based on new nodal points
                    NurbsCurve nc = NurbsCurve.Create(false, n, tempP);
                    defGeometry.Add(nc);
                }                    
                //Set output data
                DA.SetDataList(0, defGeometry);
                #endregion
            }
        }   //End of main program

        private void InterpolatePoints(Line line, int n, out List<Point3d> tempP)
        {
            tempP = new List<Point3d>(n+1);
            double[] t = LinSpace(0, 1, n + 1);
            for (int i = 0; i < t.Length; i++)
            {
                var tPm = new Point3d();
                tPm.Interpolate(line.From, line.To, t[i]);
                tempP.Add(tPm);
            }
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

        private void Shapefunctions(double L, double x, out Matrix<double> N, out Matrix<double> dN)
        {
            double N1 = -1 / L * (x - L);
            double N2 = x / L;
            double N3 = 1 - 3 * Math.Pow(x,2) / Math.Pow(L, 2) + 2 * Math.Pow(x, 3) / Math.Pow(L,3);
            double N4 = x * (1 - 2 * x / L + Math.Pow(x, 2) / Math.Pow(L, 2));
            double N5 = Math.Pow(x, 2) / Math.Pow(L, 2) * (3 - 2 * x / L);
            double N6 = Math.Pow(x, 2) / L * (x / L - 1);

            N = Matrix<double>.Build.DenseOfArray(new double[,] {
                { N1, 0, 0,  0,  0,  0, N2, 0,  0,  0,  0,  0},
                { 0, N3, 0,  0,  0, N4, 0, N5, 0,  0,  0, N6 },
                { 0, 0, N3, 0, -N4, 0, 0, 0, N5, 0, -N6, 0},
                { 0, 0, 0, N1, 0, 0, 0, 0, 0, N2, 0, 0} });

            double dN1 = -1 / L;
            double dN2 = 1 / L;
            double dN3 = -6 * x / Math.Pow(L, 2) + 6 * Math.Pow(x, 2) / Math.Pow(L,3);
            double dN4 = 1 - 4 * x / L + 3 * Math.Pow(x, 2) / Math.Pow(L, 2);
            double dN5 = 6 * x / Math.Pow(L, 2) - 6 * Math.Pow(x, 2) / Math.Pow(L,3);
            double dN6 = 3 * Math.Pow(x, 2) / Math.Pow(L, 2) - 2 * x / L;

            dN = Matrix<double>.Build.DenseOfArray(new double[,] {
            { dN1, 0, 0, 0, 0, 0, dN2, 0, 0, 0, 0, 0},
            { 0, dN3, 0, 0, 0, dN4, 0, dN5, 0, 0, 0, dN6 },
            { 0, 0, dN3, 0, -dN4, 0, 0, 0, dN5, 0, -dN6, 0},
            { 0, 0, 0, dN1, 0, 0, 0, 0, 0, dN2, 0, 0} });
        }

        private List<Point3d> CreatePointList(List<Line> geometry)
        {
            List<Point3d> points = new List<Point3d>();

            for (int i = 0; i < geometry.Count; i++) //adds every point unless it already exists in list
            {
                Line l1 = geometry[i];
                if (!points.Contains(l1.From))
                {
                    points.Add(l1.From);
                }
                if (!points.Contains(l1.To))
                {
                    points.Add(l1.To);
                }
            }

            return points;
        }

        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                return Properties.Resources.Draw;
            }
        }

        public override Guid ComponentGuid
        {
            get { return new Guid("6391b902-2ec8-487c-94fd-b921479620b3"); }
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

                Rectangle rec2 = rec1;
                rec2.X = rec1.Right + 2;

                Rectangle rec3 = rec2;
                rec3.X = rec2.Right + 2;

                Bounds = rec0;
                ButtonBounds = rec1;
                ButtonBounds2 = rec2;
                ButtonBounds3 = rec3;

            }

            GH_Palette xColor = GH_Palette.Black;
            GH_Palette yColor = GH_Palette.Black;
            GH_Palette zColor = GH_Palette.Grey;

            private Rectangle ButtonBounds { get; set; }
            private Rectangle ButtonBounds2 { get; set; }
            private Rectangle ButtonBounds3 { get; set; }
            
            protected override void Render(GH_Canvas canvas, Graphics graphics, GH_CanvasChannel channel)
            {
                base.Render(canvas, graphics, channel);
                if (channel == GH_CanvasChannel.Objects)
                {
                    GH_Capsule button = GH_Capsule.CreateTextCapsule(ButtonBounds, ButtonBounds, xColor, "Run", 3, 0);
                    button.Render(graphics, Selected, false, false);
                    button.Dispose();
                }
                if (channel == GH_CanvasChannel.Objects)
                {
                    GH_Capsule button2 = GH_Capsule.CreateTextCapsule(ButtonBounds2, ButtonBounds2, yColor, "2nd", 2, 0);
                    button2.Render(graphics, Selected, Owner.Locked, false);
                    button2.Dispose();
                }
                if (channel == GH_CanvasChannel.Objects)
                {
                    GH_Capsule button3 = GH_Capsule.CreateTextCapsule(ButtonBounds3, ButtonBounds3, zColor, "3rd", 2, 0);
                    button3.Render(graphics, Selected, Owner.Locked, false);
                    button3.Dispose();
                }
            }

            public override GH_ObjectResponse RespondToMouseDown(GH_Canvas sender, GH_CanvasMouseEvent e)
            {
                if (e.Button == MouseButtons.Left)
                {
                    RectangleF rec = ButtonBounds;
                    if (rec.Contains(e.CanvasLocation))
                    {
                        switchColor("Run");
                        if (xColor == GH_Palette.Black) { DeformedGeometry.setToggles("Run", true); }
                        if (xColor == GH_Palette.Grey) { DeformedGeometry.setToggles("Run", false); }
                        sender.Refresh();
                        return GH_ObjectResponse.Handled;
                    }
                    rec = ButtonBounds2;
                    if (rec.Contains(e.CanvasLocation))
                    {
                        if (yColor == GH_Palette.Black) { DeformedGeometry.setToggles("2nd", false); DeformedGeometry.setToggles("3rd", true); } 
                        if (yColor == GH_Palette.Grey) { DeformedGeometry.setToggles("2nd", true); DeformedGeometry.setToggles("3rd", false); }
                        switchColor("2nd");
                        switchColor("3rd");
                        sender.Refresh();
                        return GH_ObjectResponse.Handled;
                    }
                    rec = ButtonBounds3;
                    if (rec.Contains(e.CanvasLocation))
                    {
                        if (zColor == GH_Palette.Black) { DeformedGeometry.setToggles("3rd", false); DeformedGeometry.setToggles("2nd", false); }
                        if (zColor == GH_Palette.Grey) { DeformedGeometry.setToggles("3rd", true); DeformedGeometry.setToggles("2nd", true); }
                        switchColor("3rd");
                        switchColor("2nd");
                        sender.Refresh();
                        return GH_ObjectResponse.Handled;
                    }
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
                else if (button == "2nd")
                {
                    if (yColor == GH_Palette.Black) { yColor = GH_Palette.Grey; }
                    else { yColor = GH_Palette.Black; }
                }
                else if (button == "3rd")
                {
                    if (zColor == GH_Palette.Black) { zColor = GH_Palette.Grey; }
                    else { zColor = GH_Palette.Black; }
                }
            }
        }
    }
}

#region Working w/shapefunctions
//using System;
//using System.Collections.Generic;

//using Grasshopper.Kernel;
//using Rhino.Geometry;
//using System.Drawing;
//using Grasshopper.GUI.Canvas;
//using System.Windows.Forms;
//using Grasshopper.GUI;

//using MathNet.Numerics.LinearAlgebra;

//namespace Beam3D
//{
//    public class DeformedGeometry : GH_Component
//    {
//        public DeformedGeometry()
//          : base("DeformedGeometry", "DefG",
//              "Description",
//              "Koala", "3D Beam")
//        {
//        }

//        //Initialize startcondition
//        static bool startDef = true;
//        static bool p2 = true;
//        static bool p3 = false;


//        //Method to allow C# hanging of variables via GUI (see Component Visual)
//        public static void setToggles(string s, bool i)
//        {
//            if (s == "Run")
//            {
//                startDef = i;
//            }
//            else if (s == "2nd")
//            {
//                p2 = i;
//            }
//            else if (s == "3rd")
//            {
//                p3 = i;
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
//            pManager.AddNumberParameter("Deformation", "Def", "Deformations from 3DBeamCalc", GH_ParamAccess.list);
//            pManager.AddLineParameter("Geometry", "G", "Input Geometry (Line format)", GH_ParamAccess.list);
//            pManager.AddNumberParameter("Scale", "S", "The Scale Factor for Deformation", GH_ParamAccess.item, 1000);
//            pManager.AddIntegerParameter("Elements", "n", "No. of Elements", GH_ParamAccess.item, 4);
//        }

//        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
//        {
//            pManager.AddCurveParameter("Deformed Geometry", "Def.G.", "Deformed Geometry as List of Lines", GH_ParamAccess.list);
//        }

//        protected override void SolveInstance(IGH_DataAccess DA)
//        {
//            if (startDef)
//            {
//                #region Fetch input
//                //Expected inputs and outputs
//                List<double> def = new List<double>();          //input doubles of displacements/rotations
//                List<Line> geometry = new List<Line>();         //input list of original node coordinates
//                double scale = 1000;                            //input deformation scale (default == 1000)                
//                int n = 4;                                      //input no. of elements (default == 4)

//                List<Curve> defGeometry = new List<Curve>();    //output deformed geometry
//                List<Point3d> defPoints = new List<Point3d>();  //output deformed element nodes (remove?)

//                var m = Matrix<double>.Build;
//                var v = Vector<double>.Build;

//                //Set expected inputs from Indata
//                if (!DA.GetDataList(0, def)) return;
//                if (!DA.GetDataList(1, geometry)) return;
//                if (!DA.GetData(2, ref scale)) return;
//                if (!DA.GetData(3, ref n)) return;
//                #endregion

//                //List all nodes (every node only once), numbering them according to list index
//                List<Point3d> points = CreatePointList(geometry);


//                #region Create geometry (linear)
//                if (n == 1)
//                {
//                    int index = 0;
//                    //loops through all points and scales x-, y- and z-dir

//                    //u(x) = Na, N = shape func, a = nodal values (dof) 
//                    foreach (Point3d point in points)
//                    {
//                        //fetch global x,y,z placement of point
//                        double x = point.X;
//                        double y = point.Y;
//                        double z = point.Z;

//                        //scales x and z according to input Scale
//                        defPoints.Add(new Point3d(x + scale * def[index], y + scale * def[index + 1], z + scale * def[index + 2]));
//                        index += 6;
//                    }

//                    //creates deformed geometry based on initial geometry placement
//                    foreach (Line line in geometry)
//                    {
//                        //fetches index of original start and endpoint
//                        int i1 = points.IndexOf(line.From);
//                        int i2 = points.IndexOf(line.To);

//                        //creates new line based on scaled deformation of said points
//                        var tempL = new Line(defPoints[i1], defPoints[i2]);
//                        defGeometry.Add(tempL.ToNurbsCurve());
//                    }


//                    //Set output data
//                    DA.SetDataList(0, defGeometry);
//                    return; //return to grasshopper since solved linearly
//                }
//                #endregion

//                #region Create geometry
//                Matrix<double> N, B;
//                foreach (Line line in geometry)
//                {
//                    //fetches index of original start and endpoint
//                    int i1 = points.IndexOf(line.From);
//                    int i2 = points.IndexOf(line.To);

//                    //create 12x1 deformation vector for element (6dofs), scaled and populated with existing deformations
//                    var u = Vector<double>.Build.Dense(12);
//                    for (int j = 0; j < 6; j++)
//                    {
//                        u[j] = def[i1 * 6 + j];
//                        u[j + 6] = def[i2 * 6 + j];
//                    }
//                    u = scale * u;

//                    //interpolate points between line.From and line.To
//                    List<Point3d> tempP;
//                    InterpolatePoints(line, n, out tempP);


//                    double L = points[i1].DistanceTo(points[i2]);   //L is distance from startnode to endnode
//                    var x = v.Dense(n + 1);                       //maybe this should be projected x instead???
//                    for (int j = 0; j < n + 1; j++)
//                    {
//                        x[j] = j * L / n;
//                    }


//                    //Calculate 6 dofs for all new elements using shape functions (n+1 elements)
//                    Matrix<double> disp = m.Dense(n + 1, 4);
//                    Matrix<double> rot = m.Dense(n + 1, 4);

//                    for (int j = 0; j < n + 1; j++)          //x are points inbetween (?)
//                    {
//                        Shapefunctions(L, x[j], out N, out B);
//                        disp.SetRow(j, N.Multiply(u));
//                        rot.SetRow(j, B.Multiply(u));
//                    }

//                    //Calculate new nodal points
//                    for (int j = 0; j < n + 1; j++)
//                    {
//                        //original xyz                        
//                        var tP = tempP[j];

//                        //add displacement
//                        //tP.X += disp[j, 0];
//                        //tP.Y += disp[j, 1];
//                        //tP.Z += disp[j, 2];

//                        tP.X = tP.X + disp[j, 0] + tP.Z * Math.Cos(Math.PI / 2 - rot[j, 2]) + tP.Y * Math.Cos(Math.PI / 2 - rot[j, 1]);
//                        tP.Y = -Math.Cos(disp[j, 3]) * tP.Y * Math.Sin(Math.PI / 2 - rot[j, 1]) + Math.Sin(disp[j, 3]) * tP.Z + disp[j, 1];
//                        //tP.Z = disp[j, 2] + tP.Z * Math.Sin(rot[j, 2]);
//                        tP.Z = -Math.Sin(disp[j, 3]) * tP.Y + Math.Cos(disp[j, 3]) * tP.Z * Math.Sin(Math.PI / 2 - rot[j, 3]) + disp[j, 2];

//                        //replace previous xyz with displaced xyz
//                        tempP[j] = tP;
//                    }

//                    //Create Nurbscurve based on new nodal points
//                    NurbsCurve nc = NurbsCurve.Create(false, n, tempP);
//                    defGeometry.Add(nc);
//                }
//                //Set output data
//                DA.SetDataList(0, defGeometry);
//                #endregion
//            }
//        }   //End of main program

//        private void InterpolatePoints(Line line, int n, out List<Point3d> tempP)
//        {
//            tempP = new List<Point3d>(n + 1);
//            double[] t = LinSpace(0, 1, n + 1);
//            for (int i = 0; i < t.Length; i++)
//            {
//                var tPm = new Point3d();
//                tPm.Interpolate(line.From, line.To, t[i]);
//                tempP.Add(tPm);
//            }
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

//        private List<Point3d> CreatePointList(List<Line> geometry)
//        {
//            List<Point3d> points = new List<Point3d>();

//            for (int i = 0; i < geometry.Count; i++) //adds every point unless it already exists in list
//            {
//                Line l1 = geometry[i];
//                if (!points.Contains(l1.From))
//                {
//                    points.Add(l1.From);
//                }
//                if (!points.Contains(l1.To))
//                {
//                    points.Add(l1.To);
//                }
//            }

//            return points;
//        }

//        protected override System.Drawing.Bitmap Icon
//        {
//            get
//            {
//                return Properties.Resources.Draw;
//            }
//        }

//        public override Guid ComponentGuid
//        {
//            get { return new Guid("6391b902-2ec8-487c-94fd-b921479620b3"); }
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

//                Rectangle rec3 = rec2;
//                rec3.X = rec2.Right + 2;

//                Bounds = rec0;
//                ButtonBounds = rec1;
//                ButtonBounds2 = rec2;
//                ButtonBounds3 = rec3;

//            }

//            GH_Palette xColor = GH_Palette.Black;
//            GH_Palette yColor = GH_Palette.Black;
//            GH_Palette zColor = GH_Palette.Grey;

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
//                    GH_Capsule button2 = GH_Capsule.CreateTextCapsule(ButtonBounds2, ButtonBounds2, yColor, "2nd", 2, 0);
//                    button2.Render(graphics, Selected, Owner.Locked, false);
//                    button2.Dispose();
//                }
//                if (channel == GH_CanvasChannel.Objects)
//                {
//                    GH_Capsule button3 = GH_Capsule.CreateTextCapsule(ButtonBounds3, ButtonBounds3, zColor, "3rd", 2, 0);
//                    button3.Render(graphics, Selected, Owner.Locked, false);
//                    button3.Dispose();
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
//                        if (xColor == GH_Palette.Black) { DeformedGeometry.setToggles("Run", true); }
//                        if (xColor == GH_Palette.Grey) { DeformedGeometry.setToggles("Run", false); }
//                        sender.Refresh();
//                        return GH_ObjectResponse.Handled;
//                    }
//                    rec = ButtonBounds2;
//                    if (rec.Contains(e.CanvasLocation))
//                    {
//                        if (yColor == GH_Palette.Black) { DeformedGeometry.setToggles("2nd", false); DeformedGeometry.setToggles("3rd", true); }
//                        if (yColor == GH_Palette.Grey) { DeformedGeometry.setToggles("2nd", true); DeformedGeometry.setToggles("3rd", false); }
//                        switchColor("2nd");
//                        switchColor("3rd");
//                        sender.Refresh();
//                        return GH_ObjectResponse.Handled;
//                    }
//                    rec = ButtonBounds3;
//                    if (rec.Contains(e.CanvasLocation))
//                    {
//                        if (zColor == GH_Palette.Black) { DeformedGeometry.setToggles("3rd", false); DeformedGeometry.setToggles("2nd", false); }
//                        if (zColor == GH_Palette.Grey) { DeformedGeometry.setToggles("3rd", true); DeformedGeometry.setToggles("2nd", true); }
//                        switchColor("3rd");
//                        switchColor("2nd");
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
//                else if (button == "2nd")
//                {
//                    if (yColor == GH_Palette.Black) { yColor = GH_Palette.Grey; }
//                    else { yColor = GH_Palette.Black; }
//                }
//                else if (button == "3rd")
//                {
//                    if (zColor == GH_Palette.Black) { zColor = GH_Palette.Grey; }
//                    else { zColor = GH_Palette.Black; }
//                }
//            }
//        }
//    }
//}
#endregion