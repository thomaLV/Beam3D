using System;
using System.Collections.Generic;

using Grasshopper.Kernel;
using Rhino.Geometry;
using System.Drawing;
using Grasshopper.GUI.Canvas;
using System.Windows.Forms;
using Grasshopper.GUI;

namespace Beam3D
{
    public class SetMoments : GH_Component
    {
        public SetMoments()
          : base("SetMoments", "Nickname",
                  "Description",
                  "Koala", "3D Beam")
        {
        }
        //Initialize moments
        static int mx = 1;
        static int my = 1;
        static int mz = 1;


        //Method to allow c hanging of variables via GUI (see Component Visual)
        public static void setMom(string s, int i)
        {
            if (s == "MX")
            {
                mx = i;
            }
            else if (s == "MY")
            {
                my = i;
            }
            else if (s == "MZ")
            {
                mz = i;
            }
        }

        public override void CreateAttributes()
        {
            m_attributes = new Attributes_Custom(this);
        }

        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddPointParameter("Points", "P", "Points to apply moment", GH_ParamAccess.list);
            pManager.AddNumberParameter("Moment", "M", "Moment Magnitude [kNm]", GH_ParamAccess.list);
        }

        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddTextParameter("MomentLoads", "ML", "MomentLoads formatted for Beam Calculation", GH_ParamAccess.list);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            #region Fetch inputs
            //Expected inputs and output
            List<Point3d> pointList = new List<Point3d>();              //List of points where load will be applied
            List<double> momentList = new List<double>();                 //List or value of load applied
            List<string> pointInStringFormat = new List<string>();      //preallocate final string output

            //Set expected inputs from Indata
            if (!DA.GetDataList(0, pointList)) return;
            if (!DA.GetDataList(1, momentList)) return;
            #endregion

            #region Format output
            string vectorString;

            for (int i = 0, j = 0; i < pointList.Count; i++)   //Format stringline for all points (identical boundary conditions for all points)
            {
                    vectorString = momentList[j] * mx + "," + momentList[j] * my + "," + momentList[j] * mz;
                    pointInStringFormat.Add(pointList[i].X + "," + pointList[i].Y + "," + pointList[i].Z + ":" + vectorString);
                if (j < momentList.Count - 1)
                {
                    j++;
                }
            }
            #endregion

            //Set output data
            DA.SetDataList(0, pointInStringFormat);
        }

        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                return Properties.Resources.Moments;
            }
        }

        public override Guid ComponentGuid
        {
            get { return new Guid("540c5cd8-b017-45d3-b3d1-cb1bf0c9051c"); }
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

                BoundsAllButtons = rec0;
                Bounds = rec0;
                ButtonBounds = rec1;
                ButtonBounds2 = rec2;
                ButtonBounds3 = rec3;
            }

            GH_Palette xColor = GH_Palette.Black;
            GH_Palette yColor = GH_Palette.Black;
            GH_Palette zColor = GH_Palette.Black;

            private Rectangle BoundsAllButtons { get; set; }
            private Rectangle ButtonBounds { get; set; }
            private Rectangle ButtonBounds2 { get; set; }
            private Rectangle ButtonBounds3 { get; set; }

            protected override void Render(GH_Canvas canvas, Graphics graphics, GH_CanvasChannel channel)
            {
                base.Render(canvas, graphics, channel);
                if (channel == GH_CanvasChannel.Objects)
                {
                    GH_Capsule button = GH_Capsule.CreateTextCapsule(ButtonBounds, ButtonBounds, xColor, "MX", 2, 0);
                    button.Render(graphics, Selected, Owner.Locked, false);
                    button.Dispose();
                }
                if (channel == GH_CanvasChannel.Objects)
                {
                    GH_Capsule button2 = GH_Capsule.CreateTextCapsule(ButtonBounds2, ButtonBounds2, yColor, "MY", 2, 0);
                    button2.Render(graphics, Selected, Owner.Locked, false);
                    button2.Dispose();
                }
                if (channel == GH_CanvasChannel.Objects)
                {
                    GH_Capsule button3 = GH_Capsule.CreateTextCapsule(ButtonBounds3, ButtonBounds3, zColor, "MZ", 2, 0);
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
                        switchColor("MX");
                    }
                    rec = ButtonBounds2;
                    if (rec.Contains(e.CanvasLocation))
                    {
                        switchColor("MY");
                    }
                    rec = ButtonBounds3;
                    if (rec.Contains(e.CanvasLocation))
                    {
                        switchColor("MZ");
                    }
                    rec = BoundsAllButtons;
                    if (rec.Contains(e.CanvasLocation))
                    {
                        if (xColor == GH_Palette.Black) { setMom("MX", 0); }
                        if (xColor == GH_Palette.Grey) { setMom("MX", 1); }
                        if (yColor == GH_Palette.Black) { setMom("MY", 0); }
                        if (yColor == GH_Palette.Grey) { setMom("MY", 1); }
                        if (zColor == GH_Palette.Black) { setMom("MZ", 0); }
                        if (zColor == GH_Palette.Grey) { setMom("MZ", 1); }
                        Owner.ExpireSolution(true);
                    }
                    return GH_ObjectResponse.Handled;
                }
                return base.RespondToMouseDown(sender, e);
            }

            private void switchColor(string button)
            {
                if (button == "MX")
                {
                    if (xColor == GH_Palette.Black) { xColor = GH_Palette.Grey; }
                    else { xColor = GH_Palette.Black; }
                }
                else if (button == "MY")
                {
                    if (yColor == GH_Palette.Black) { yColor = GH_Palette.Grey; }
                    else { yColor = GH_Palette.Black; }
                }
                else if (button == "MZ")
                {
                    if (yColor == GH_Palette.Black) { yColor = GH_Palette.Grey; }
                    else { yColor = GH_Palette.Black; }
                }
                Owner.ExpireSolution(true);
            }
        }
    }
}
