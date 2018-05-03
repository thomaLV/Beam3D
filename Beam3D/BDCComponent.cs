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
    public class BDCComponent : GH_Component
    {
        public BDCComponent()
          : base("BDCComponent", "BDCs",
              "Description",
              "Koala", "3D Beam")
        {
        }

        //Initialize BDCs
        private static int x = 0;
        private static int y = 0;
        private static int z = 0;
        private static int rx = 0;
        private static int ry = 0;
        private static int rz = 0;

        public static void setBDC(string s, int i)
        {
            if (s == "X")
            {
                x = i;
            }
            else if (s == "Y")
            {
                y = i;
            }
            else if (s == "Z")
            {
                z = i;
            }
            else if (s == "RX")
            {
                rx = i;
            }
            else if (s == "RY")
            {
                ry = i;
            }
            else if (s == "RZ")
            {
                rz = i;
            }
            //Grasshopper.Instances.ActiveCanvas.Document.ExpireSolution();
            //Grasshopper.Instances.ActiveCanvas.Document.NewSolution(false);
        }

        public override void CreateAttributes()
        {
            m_attributes = new Attributes_Custom(this);
        }

        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddPointParameter("Points", "P", "Points to apply Boundary Conditions", GH_ParamAccess.list);
            //pManager.AddLineParameter("Geometry", "G", "Geometry", GH_ParamAccess.list);
            //pManager.AddIntegerParameter("Boundary Conditions", "BDC", "Boundary Conditions translation of x,y,z,rx,ry,rz where 0=clamped and 1=free", GH_ParamAccess.list, new List<int>());
            //pManager.AddTextParameter("Locked direction", "Ldir", "Lock x, y or z direction for deformation", GH_ParamAccess.item, "");
        }

        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddTextParameter("B.Cond.", "BDC", "Boundary Conditions for 3D Beam Calculation", GH_ParamAccess.list);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            #region Fetch inputs
            //Expected inputs
            List<Point3d> pointList = new List<Point3d>();          //List of points where BDC is to be applied
            List<string> pointInStringFormat = new List<string>();  //output in form of list of strings


            //Set expected inputs from Indata and aborts with error message if input is incorrect
            if (!DA.GetDataList(0, pointList)) return;
            #endregion
            
            #region Format output
            string BDCString = x + "," + y + "," + z + "," + rx + "," + ry + "," + rz;

            for (int i = 0; i < pointList.Count; i++)   //Format stringline for all points (identical boundary conditions for all points)
            {
                pointInStringFormat.Add(pointList[i].X + "," + pointList[i].Y + "," + pointList[i].Z + ":" + BDCString);
            }
            #endregion

            #region TODO (implement lock_dir)
            //lock_dir = lock_dir.ToUpper();

            //Preallocate temporary variables
            //string BDCString;
            //int tbdcx = 0;
            //int tbdcy = 0;
            //int tbdcz = 0;
            //int rbdcx = 0;
            //int rbdcy = 0;
            //int rbdcz = 0;


            //if (lock_dir == "")
            //{
            //    if (BDC.Count == 1) //Boundary condition input for identical conditions in all points. Split into if/else for optimization
            //    {
            //        tbdcx = BDC[0];
            //        tbdcy = BDC[0];
            //        tbdcz = BDC[0];
            //        rbdcx = BDC[0];
            //        rbdcy = BDC[0];
            //        rbdcz = BDC[0];

            //        BDCString = tbdcx + "," + tbdcy + "," + tbdcz + "," + rbdcx + "," + rbdcy + "," + rbdcz;
            //        for (int i = 0; i < pointList.Count; i++)   //Format stringline for all points (identical boundary conditions for all points)
            //        {
            //            pointInStringFormat.Add(pointList[i].X + "," + pointList[i].Y + "," + pointList[i].Z + ":" + BDCString);
            //        }
            //    }
            //    else if (BDC.Count == 6) //Boundary condition input for identical conditions in all points. Split into if/else for optimization
            //    {
            //        tbdcx = BDC[0];
            //        tbdcy = BDC[1];
            //        tbdcz = BDC[2];
            //        rbdcx = BDC[3];
            //        rbdcy = BDC[4];
            //        rbdcz = BDC[5];

            //        BDCString = tbdcx + "," + tbdcy + "," + tbdcz + "," + rbdcx + "," + rbdcy + "," + rbdcz;

            //        for (int i = 0; i < pointList.Count; i++)   //Format stringline for all points (identical boundary conditions for all points)
            //        {
            //            pointInStringFormat.Add(pointList[i].X + "," + pointList[i].Y + "," + pointList[i].Z + ":" + BDCString);
            //        }
            //    }
            //    else    //BDCs are not identical for all points
            //    {
            //        for (int i = 0; i < pointList.Count; i++)
            //        {
            //            if (i > (BDC.Count / 6) - 1)  //Are there more points than BDCs given? (BDC always lists x,y,z per point)
            //            {
            //                //use values from last BDC in list of BDCs
            //                BDCString = tbdcx + "," + tbdcy + "," + tbdcz + "," + rbdcx + "," + rbdcy + "," + rbdcz;
            //            }
            //            else
            //            {
            //                //retrieve BDC for x,y,z-dir
            //                tbdcx = BDC[i * 6 + 0];
            //                tbdcy = BDC[i * 6 + 1];
            //                tbdcz = BDC[i * 6 + 2];
            //                rbdcx = BDC[i * 6 + 3];
            //                rbdcy = BDC[i * 6 + 4];
            //                rbdcz = BDC[i * 6 + 5];
            //                BDCString = tbdcx + "," + tbdcy + "," + tbdcz + "," + rbdcx + "," + rbdcy + "," + rbdcz;
            //            }
            //            pointInStringFormat.Add(pointList[i].X + "," + pointList[i].Y + "," + pointList[i].Z + ":" + BDCString);    //Add stringline to list of strings
            //        }
            //    }
            //}
            //else
            //{
            //    bool lx = false;
            //    bool ly = false;
            //    bool lz = false;
            //    bool rx = false;
            //    bool ry = false;
            //    bool rz = false;

            //    if (lock_dir == "X")
            //    {
            //        lx = true;
            //        tbdcx = 0;
            //    }
            //    else if (lock_dir == "Y")
            //    {
            //        ly = true;
            //        tbdcy = 0;
            //    }
            //    else if (lock_dir == "Z")
            //    {
            //        ly = true;
            //        tbdcz = 0;
            //    }
            //    else if (lock_dir == "MX")
            //    {
            //        rx = true;
            //        rbdcx = 0;
            //    }
            //    else if (lock_dir == "MY")
            //    {
            //        ry = true;
            //        rbdcy = 0;
            //    }
            //    else if (lock_dir == "MZ")
            //    {
            //        ry = true;
            //        rbdcz = 0;
            //    }

            //    List<Point3d> points = CreatePointList(geometry);
            //    for (int i = 0; i < pointList.Count; i++)
            //    {
            //        points.Remove(pointList[i]);
            //    }

            //    for (int i = 0; i < points.Count; i++)
            //    {
            //        if (!lx) tbdcx = 1;
            //        if (!ly) tbdcy = 1;
            //        if (!lz) tbdcz = 1;
            //        if (!rx) rbdcx = 1;
            //        if (!ry) rbdcy = 1;
            //        if (!rz) rbdcz = 1;

            //        BDCString = tbdcx + "," + tbdcy + "," + tbdcz + "," + rbdcx + "," + rbdcy + "," + rbdcz;
            //        pointInStringFormat.Add(points[i].X + "," + points[i].Y + "," + points[i].Z + ":" + BDCString);
            //    }

            //    if (BDC.Count == 1) //Boundary condition input for identical conditions in all points. Split into if/else for optimization
            //    {
            //        if (!lx) tbdcx = BDC[0];
            //        if (!ly) tbdcy = BDC[0];
            //        if (!lz) tbdcz = BDC[0];
            //        if (!rx) rbdcx = BDC[0];
            //        if (!ry) rbdcy = BDC[0];
            //        if (!rz) rbdcz = BDC[0];

            //        BDCString = tbdcx + "," + tbdcy + "," + tbdcz + "," + rbdcx + "," + rbdcy + "," + rbdcz;
            //        for (int i = 0; i < pointList.Count; i++)   //Format stringline for all points (identical boundary conditions for all points)
            //        {
            //            pointInStringFormat.Add(pointList[i].X + "," + pointList[i].Y + "," + pointList[i].Z + ":" + BDCString);
            //        }
            //    }
            //    else if (BDC.Count == 6) //Boundary condition input for identical conditions in all points. Split into if/else for optimization
            //    {
            //        if (!lx) tbdcx = BDC[0];
            //        if (!ly) tbdcy = BDC[1];
            //        if (!lz) tbdcz = BDC[2];
            //        if (!rx) rbdcx = BDC[3];
            //        if (!ry) rbdcy = BDC[4];
            //        if (!rz) rbdcz = BDC[5];
            //        BDCString = tbdcx + "," + tbdcy + "," + tbdcz + "," + rbdcx + "," + rbdcy + "," + rbdcz;
            //        for (int i = 0; i < pointList.Count; i++)   //Format stringline for all points (identical boundary conditions for all points)
            //        {
            //            pointInStringFormat.Add(pointList[i].X + "," + pointList[i].Y + "," + pointList[i].Z + ":" + BDCString);
            //        }
            //    }
            //    else    //BDCs are not identical for all points
            //    {
            //        for (int i = 0; i < pointList.Count; i++)
            //        {
            //            if (i > (BDC.Count / 6) - 1)  //Are there more points than BDCs given? (BDC always lists x,y,z per point)
            //            {
            //                BDCString = tbdcx + "," + tbdcy + "," + tbdcz + "," + rbdcx + "," + rbdcy + "," + rbdcz;
            //            }
            //            else
            //            {
            //                //retrieve BDC for x,y,z-dir
            //                if (!lx) tbdcx = BDC[i * 6];
            //                if (!ly) tbdcy = BDC[i * 6 + 1];
            //                if (!lz) tbdcz = BDC[i * 6 + 2];
            //                if (!rx) rbdcx = BDC[i * 6 + 3];
            //                if (!ry) rbdcy = BDC[i * 6 + 4];
            //                if (!rz) rbdcz = BDC[i * 6 + 5];
            //                BDCString = tbdcx + "," + tbdcy + "," + tbdcz + "," + rbdcx + "," + rbdcy + "," + rbdcz;
            //            }
            //            pointInStringFormat.Add(pointList[i].X + "," + pointList[i].Y + "," + pointList[i].Z + ":" + BDCString);    //Add stringline to list of strings
            //        }
            //    }
            //}
            #endregion

            DA.SetDataList(0, pointInStringFormat);
        } //End of main program

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

        //TODO (implement lock_dir buttons)
        //public static void switchBDC(string button)
        //{
        //    if (button == "X")
        //    {
        //        if (BoundaryConditions3DBeamComponent.) { xColor = GH_Palette.Black; }
        //    }
        //    else if (button == "Y")
        //    {
        //        if (yColor == GH_Palette.Black) { yColor = GH_Palette.Grey; }
        //        else { yColor = GH_Palette.Black; }
        //    }
        //    else if (button == "Z")
        //    {
        //        if (zColor == GH_Palette.Black) { zColor = GH_Palette.Grey; }
        //        else { zColor = GH_Palette.Black; }
        //    }
        //    else if (button == "RX")
        //    {
        //        if (rxColor == GH_Palette.Black) { rxColor = GH_Palette.Grey; }
        //        else { rxColor = GH_Palette.Black; }
        //    }
        //    else if (button == "RY")
        //    {
        //        if (ryColor == GH_Palette.Black) { ryColor = GH_Palette.Grey; }
        //        else { ryColor = GH_Palette.Black; }
        //    }
        //    else if (button == "RZ")
        //    {
        //        if (rzColor == GH_Palette.Black) { rzColor = GH_Palette.Grey; }
        //        else { rzColor = GH_Palette.Black; }
        //    }
        //}

        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                return Properties.Resources.BDCs;
            }
        }
        
        public override Guid ComponentGuid
        {
            get { return new Guid("c9c208e0-b10b-4ecb-a5ef-57d86a4df109"); }
        }
    }


    /// Component Visual//
    public class Attributes_Custom : Grasshopper.Kernel.Attributes.GH_ComponentAttributes
    {
        public Attributes_Custom(GH_Component owner) : base(owner) { }
        protected override void Layout()
        {
            base.Layout();

            Rectangle rec0 = GH_Convert.ToRectangle(Bounds);

            rec0.Height += 42;

            Rectangle rec1 = rec0;
            rec1.X = rec0.Left + 1;
            rec1.Y = rec0.Bottom - 42;
            rec1.Width = (rec0.Width) / 3 + 1;
            rec1.Height = 22;
            rec1.Inflate(-2, -2);

            Rectangle rec2 = rec1;
            rec2.X = rec1.Right + 2;

            Rectangle rec3 = rec2;
            rec3.X = rec2.Right + 2;

            Rectangle rec4 = rec1;
            rec4.Y = rec1.Bottom + 2;

            Rectangle rec5 = rec4;
            rec5.X = rec4.Right + 2;

            Rectangle rec6 = rec5;
            rec6.X = rec2.Right + 2;

            Bounds = rec0;
            ButtonBounds = rec1;
            ButtonBounds2 = rec2;
            ButtonBounds3 = rec3;
            ButtonBounds4 = rec4;
            ButtonBounds5 = rec5;
            ButtonBounds6 = rec6;

        }

        GH_Palette xColor = GH_Palette.Black;
        GH_Palette yColor = GH_Palette.Black;
        GH_Palette zColor = GH_Palette.Black;
        GH_Palette rxColor = GH_Palette.Black;
        GH_Palette ryColor = GH_Palette.Black;
        GH_Palette rzColor = GH_Palette.Black;

        private Rectangle ButtonBounds { get; set; }
        private Rectangle ButtonBounds2 { get; set; }
        private Rectangle ButtonBounds3 { get; set; }
        private Rectangle ButtonBounds4 { get; set; }
        private Rectangle ButtonBounds5 { get; set; }
        private Rectangle ButtonBounds6 { get; set; }

        protected override void Render(GH_Canvas canvas, Graphics graphics, GH_CanvasChannel channel)
        {
            base.Render(canvas, graphics, channel);
            if (channel == GH_CanvasChannel.Objects)
            {
                GH_Capsule button = GH_Capsule.CreateTextCapsule(ButtonBounds, ButtonBounds, xColor, "X", 3, 0);
                button.Render(graphics, Selected, false, false);
                button.Dispose();
            }
            if (channel == GH_CanvasChannel.Objects)
            {
                GH_Capsule button2 = GH_Capsule.CreateTextCapsule(ButtonBounds2, ButtonBounds2, yColor, "Y", 2, 0);
                button2.Render(graphics, Selected, Owner.Locked, false);
                button2.Dispose();
            }
            if (channel == GH_CanvasChannel.Objects)
            {
                GH_Capsule button3 = GH_Capsule.CreateTextCapsule(ButtonBounds3, ButtonBounds3, zColor, "Z", 2, 0);
                button3.Render(graphics, Selected, Owner.Locked, false);
                button3.Dispose();
            }
            if (channel == GH_CanvasChannel.Objects)
            {
                GH_Capsule button4 = GH_Capsule.CreateTextCapsule(ButtonBounds4, ButtonBounds4, rxColor, "RX", 2, 0);
                button4.Render(graphics, Selected, Owner.Locked, false);
                button4.Dispose();
            }
            if (channel == GH_CanvasChannel.Objects)
            {
                GH_Capsule button5 = GH_Capsule.CreateTextCapsule(ButtonBounds5, ButtonBounds5, ryColor, "RY", 2, 0);
                button5.Render(graphics, Selected, Owner.Locked, false);
                button5.Dispose();
            }
            if (channel == GH_CanvasChannel.Objects)
            {
                GH_Capsule button6 = GH_Capsule.CreateTextCapsule(ButtonBounds6, ButtonBounds6, rzColor, "RZ", 2, 0);
                button6.Render(graphics, Selected, Owner.Locked, false);
                button6.Dispose();
            }
        }

        public override GH_ObjectResponse RespondToMouseDown(GH_Canvas sender, GH_CanvasMouseEvent e)
        {
            if (e.Button == MouseButtons.Left)
            {
                RectangleF rec = ButtonBounds;
                if (rec.Contains(e.CanvasLocation))
                {
                    switchColor("X");
                    if (xColor == GH_Palette.Black) { BDCComponent.setBDC("X", 0); }
                    if (xColor == GH_Palette.Grey) { BDCComponent.setBDC("X", 1); }
                    sender.Refresh();
                    //return GH_ObjectResponse.Handled;
                }
                rec = ButtonBounds2;
                if (rec.Contains(e.CanvasLocation))
                {
                    switchColor("Y");
                    if (yColor == GH_Palette.Black) { BDCComponent.setBDC("Y", 0); }
                    if (yColor == GH_Palette.Grey) { BDCComponent.setBDC("Y", 1); }
                    sender.Refresh();
                    //return GH_ObjectResponse.Handled;
                }
                rec = ButtonBounds3;
                if (rec.Contains(e.CanvasLocation))
                {
                    switchColor("Z");
                    if (zColor == GH_Palette.Black) { BDCComponent.setBDC("Z", 0); }
                    if (zColor == GH_Palette.Grey) { BDCComponent.setBDC("Z", 1); }
                    sender.Refresh();
                    //return GH_ObjectResponse.Handled;
                }
                rec = ButtonBounds4;
                if (rec.Contains(e.CanvasLocation))
                {
                    switchColor("RX");
                    if (rxColor == GH_Palette.Black) { BDCComponent.setBDC("RX", 0); }
                    if (rxColor == GH_Palette.Grey) { BDCComponent.setBDC("RX", 1); }
                    sender.Refresh();
                    //return GH_ObjectResponse.Handled;
                }
                rec = ButtonBounds5;
                if (rec.Contains(e.CanvasLocation))
                {
                    switchColor("RY");
                    if (ryColor == GH_Palette.Black) { BDCComponent.setBDC("RY", 0); }
                    if (ryColor == GH_Palette.Grey) { BDCComponent.setBDC("RY", 1); }
                    sender.Refresh();
                    //return GH_ObjectResponse.Handled;
                }
                rec = ButtonBounds6;
                if (rec.Contains(e.CanvasLocation))
                {
                    switchColor("RZ");
                    if (rzColor == GH_Palette.Black) { BDCComponent.setBDC("RZ", 0); }
                    if (rzColor == GH_Palette.Grey) { BDCComponent.setBDC("RZ", 1); }
                    sender.Refresh();
                    //return GH_ObjectResponse.Handled;
                }
            }
            Owner.ExpireSolution(true);
            return base.RespondToMouseDown(sender, e);
        }

        private void switchColor(string button)
        {
            if (button == "X")
            {
                if (xColor == GH_Palette.Black) { xColor = GH_Palette.Grey; }
                else { xColor = GH_Palette.Black; }
            }
            else if (button == "Y")
            {
                if (yColor == GH_Palette.Black) { yColor = GH_Palette.Grey; }
                else { yColor = GH_Palette.Black; }
            }
            else if (button == "Z")
            {
                if (zColor == GH_Palette.Black) { zColor = GH_Palette.Grey; }
                else { zColor = GH_Palette.Black; }
            }
            else if (button == "RX")
            {
                if (rxColor == GH_Palette.Black) { rxColor = GH_Palette.Grey; }
                else { rxColor = GH_Palette.Black; }
            }
            else if (button == "RY")
            {
                if (ryColor == GH_Palette.Black) { ryColor = GH_Palette.Grey; }
                else { ryColor = GH_Palette.Black; }
            }
            else if (button == "RZ")
            {
                if (rzColor == GH_Palette.Black) { rzColor = GH_Palette.Grey; }
                else { rzColor = GH_Palette.Black; }
            }
        }
    }

}