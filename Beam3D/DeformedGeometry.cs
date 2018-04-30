﻿using System;
using System.Collections.Generic;

using Grasshopper.Kernel;
using Rhino.Geometry;
using System.Drawing;
using Grasshopper.GUI.Canvas;
using System.Windows.Forms;
using Grasshopper.GUI;

using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;

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

                //number of divisions(?)

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

                ////polynomial order (should be 1, 2??)
                //int p;
                //if (p2)
                //{
                //     p = 2;
                //}
                //else
                //{
                //    p = 3;
                //}

                if (n == 1)
                {
                    //Set output data
                    DA.SetDataList(0, defGeometry);
                }

                Matrix<double> N, dN;


                #region Calculate deformed geometry
                foreach (Line line in geometry)
                {
                    //fetches index of original start and endpoint
                    int i1 = points.IndexOf(line.From);
                    int i2 = points.IndexOf(line.To);

                    //create 12x1 deformation vector for element (6dofs), scaled and populated with existing deformations
                    var u = Vector<double>.Build.Dense(12);
                    for (int j = 0; j < 6; j++)
                    {
                        u[j    ] = def[i1 + j];
                        u[j + 6] = def[i2 + j];
                    }
                    u = scale * u;

                    //interpolate points between line.From and line.To
                    List<Point3d> tempP = new List<Point3d>(n + 1);
                    var tPm = new Point3d(); var tPse = new Point3d(); var tPla = new Point3d();
                    tPm.Interpolate(line.From, line.To, 0.5);
                    tempP.Add(line.From);
                    tPse.Interpolate(line.From, tPm, 0.5);
                    tPla.Interpolate(tPm, line.To, 0.5);
                    tempP.Add(tPse);
                    tempP.Add(tPm);
                    tempP.Add(tPla);
                    tempP.Add(line.To);


                    double L = points[i1].DistanceTo(points[i2]);   //L is distance from startnode to endnode
                    var x = v.Dense(n + 1);                       //maybe this should be projected x instead???
                    for (int j = 0; j < n + 1; j++)
                    {
                        x[j] = j * L / n;
                    }


                    //Calculate 6 dofs for all new elements using shape functions (n+1 elements)
                    Matrix<double> disp = m.Dense(n+1, 4);
                    Matrix<double> rot = m.Dense(n+1, 4);
                    for (int j = 0; j < n+1; j++)          //x are points inbetween (?)
                    {
                        Shapefunctions(L, x[j], out N, out dN);
                        disp.SetRow(j, N.Multiply(u));
//                        rot.SetRow(j, dN.Multiply(u));
                    }

                    //Calculate new nodal points
                    for (int j = 0; j < n+1; j++)
                    {
                        //original xyz                        
                        var tP = tempP[j];

                        //add displacement
                        tP.X += disp[j, 1];
                        tP.Y += disp[j, 2];
                        tP.Z += disp[j, 3];

                        //replace previous xyz with displaced xyz
                        tempP[j] = tP;
                    }

                    //Create Nurbscurve based on new nodal points
                    NurbsCurve nc = NurbsCurve.Create(false, 3, tempP);
                    defGeometry.Add(nc);
                }
                #endregion
                

                //Set output data
                DA.SetDataList(0, defGeometry);
            }
        }   //End of main program


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