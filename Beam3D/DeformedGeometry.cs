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

        ////Initialize startcondition
        //static bool startDef = true;


        ////Method to allow C# hanging of variables via GUI (see Component Visual)
        //public static void setToggles(string s, bool i)
        //{
        //    if (s == "Color")
        //    {
        //        startDef = i;
        //    }
        //}

        //public override void CreateAttributes()
        //{
        //    m_attributes = new Attributes_Custom(this);
        //}

        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddNumberParameter("Stress", "Ss", "Nodal stress", GH_ParamAccess.list);
            pManager.AddNumberParameter("Strain", "Sn", "Nodal strain", GH_ParamAccess.list);
            pManager.AddGenericParameter("Matrix Deformations", "MDef", "Matrix Deformations from 3DBeamCalc", GH_ParamAccess.item);
            pManager.AddPointParameter("New base points", "NBP", "New base points from Calc component", GH_ParamAccess.list);
            pManager.AddNumberParameter("Scale", "S", "The Scale Factor for Deformation", GH_ParamAccess.item, 1000);
        }

        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddNumberParameter("Pure Axial stress", "PA SS", "Pure axial stress per sub-element", GH_ParamAccess.list);
            pManager.AddNumberParameter("Pure Axial strain", "PA SN", "Pure axial strain per sub-element", GH_ParamAccess.list);
            pManager.AddNumberParameter("Axial stress", "A SS", "Axial stress per sub-element", GH_ParamAccess.list);
            pManager.AddNumberParameter("Axial strain", "A SN", "Axial strain per sub-element", GH_ParamAccess.list);
            pManager.AddCurveParameter("Deformed Geometry", "Def.G.", "Deformed Geometry as List of Lines", GH_ParamAccess.list);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            #region Fetch
            //Expected inputs and outputs
            List<Curve> defC = new List<Curve>();
            List<double> stress = new List<double>();
            List<double> strain = new List<double>();
            Matrix<double> def = Matrix<double>.Build.Dense(1, 1);
            List<Point3d> oldXYZ = new List<Point3d>();
            double scale = 1000; //input deformation scale


            //Set expected inputs from Indata
            if (!DA.GetDataList(0, stress)) return;
            if (!DA.GetDataList(1, strain)) return;
            if (!DA.GetData(2, ref def)) return;
            if (!DA.GetDataList(3, oldXYZ)) return;
            if (!DA.GetData(4, ref scale)) return;
            #endregion

            #region Deformed geometry
            //no. of nodes per element
            int n = def.ColumnCount / 6;
            int ns = n - 1;

            //scale deformations
            def = scale * def;

            if (oldXYZ.Count == 0) return;
            //Calculate new nodal points
            for (int i = 0; i < def.RowCount; i++)
            {
                List<Point3d> tempNew = new List<Point3d>();
                for (int j = 0; j < n; j++)
                {
                    //original xyz
                    var tP = oldXYZ[i * n + j];

                    //add deformations
                    tP.X = tP.X + def[i, j * 6];
                    tP.Y = tP.Y + def[i, j * 6 + 1];
                    tP.Z = tP.Z + def[i, j * 6 + 2];

                    //replace previous xyz with displaced xyz
                    tempNew.Add(tP);
                }
                //Create Curve based on new nodal points(degree = 3)
                Curve nc = Curve.CreateInterpolatedCurve(tempNew, 3);
                defC.Add(nc);
            }
            #endregion

            List<double> ss_x = new List<double>();
            List<double> sn_x = new List<double>();
            List<double> ss_y = new List<double>();
            List<double> sn_y = new List<double>();
            List<double> ss_z = new List<double>();
            List<double> sn_z = new List<double>();

            for (int i = 0; i < stress.Count / 3; i++)
            {
                ss_x.Add(stress[i * 3]);
                sn_x.Add(strain[i * 3]);
                ss_y.Add(stress[i * 3 + 1]);
                sn_y.Add(strain[i * 3 + 1]);
                ss_z.Add(stress[i * 3 + 2]);
                sn_z.Add(strain[i * 3 + 2]);
            }

            ss_x = GetAverage(ss_x, ns, defC.Count);
            sn_x = GetAverage(sn_x, ns, defC.Count);
            ss_y = GetAverage(ss_y, ns, defC.Count);
            sn_y = GetAverage(sn_y, ns, defC.Count);
            ss_z = GetAverage(ss_z, ns, defC.Count);
            sn_z = GetAverage(sn_z, ns, defC.Count);

            List<double> ss = new List<double>();
            List<double> sn = new List<double>();

            for (int i = 0; i < ss_x.Count; i++)
            {
                if (ss_x[i] > 0)
                {
                    ss.Add(ss_x[i] + Math.Abs(ss_y[i]) + Math.Abs(ss_z[i]));
                    sn.Add(ss_x[i] + Math.Abs(sn_y[i]) + Math.Abs(sn_z[i]));
                }
                else
                {
                    ss.Add(ss_x[i] - Math.Abs(ss_y[i]) - Math.Abs(ss_z[i]));
                    sn.Add(sn_x[i] - Math.Abs(sn_y[i]) - Math.Abs(sn_z[i]));
                }
            }



            DA.SetDataList(0, ss_x);
            DA.SetDataList(1, sn_x);
            DA.SetDataList(2, ss);
            DA.SetDataList(3, sn);
            DA.SetDataList(4, defC);
        }//End of main program

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

        private List<double> GetAverage(List<double> s, int n, int el)
        {
            var s_avg = new List<double>();
            for (int i = 0, ct = 0; s_avg.Count < el*n; i++)
            {
                if (ct == n)
                {
                    ct = 0;
                    continue;
                }
                s_avg.Add((s[i] + s[i + 1]) / 2);
                ct++;
            }
            return s_avg;
        }
    }
}