using System;
using System.Collections.Generic;

using Grasshopper.Kernel;
using Rhino.Geometry;

namespace Beam3D
{
    public class Avg_element : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the Avg_element class.
        /// </summary>
        public Avg_element()
          : base("Avg_sub_element", "AvgSE",
              "Description",
              "Koala", "3D Beam")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddPointParameter("New base nodes", "NBN", "New base nodes", GH_ParamAccess.list);
            pManager.AddCurveParameter("Deformed geometry", "DG", "Deformed geometry", GH_ParamAccess.list);
            pManager.AddNumberParameter("Stress", "Ss", "Nodal stress", GH_ParamAccess.list);
            pManager.AddNumberParameter("Strain", "Sn", "Nodal strain", GH_ParamAccess.list);
            pManager.AddIntegerParameter("Elements", "E", "Number of sub-elements", GH_ParamAccess.item, 1);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddNumberParameter("Stress", "A", "Stress per sub-element", GH_ParamAccess.list);
            pManager.AddNumberParameter("Strain", "A", "Strain per sub-element", GH_ParamAccess.list);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            #region Fetch inputs
            //Expected inputs
            List<Point3d> pointList = new List<Point3d>();          //List of points where BDC is to be applied
            List<Curve> curveList = new List<Curve>();  //output in form of list of strings
            List<double> stress = new List<double>();  //output in form of list of strings
            List<double> strain = new List<double>();  //output in form of list of strings
            int n = 1;  //output in form of list of strings

            List<double> ss = new List<double>();
            List<double> sn = new List<double>();

            //Set expected inputs from Indata and aborts with error message if input is incorrect
            if (!DA.GetDataList(0, pointList)) return;
            if (!DA.GetDataList(1, curveList)) return;
            if (!DA.GetDataList(2, stress)) return;
            if (!DA.GetDataList(3, strain)) return;
            if (!DA.GetData(4, ref n)) return;
            #endregion

            List<double> ss_x = new List<double>();
            List<double> sn_x = new List<double>();
            
            for (int str = 0; str < 1; str++)
            {
                List<double> sst = new List<double>();
                List<double> snt = new List<double>();
                for (int i = str; i < stress.Count/3; i++)
                {
                    sst.Add(stress[i * 3]);
                    snt.Add(strain[i * 3]);
                }

                ss_x = GetAverage(sst, n);
                sn_x = GetAverage(snt, n);
            }
            


            DA.SetDataList(0, ss_x);
            DA.SetDataList(1, sn_x);
        }

        private List<double> GetAverage(List<double> s, int n)
        {
            var s_avg = new List<double>();
            for (int i = 0, ct = 0; s_avg.Count < s.Count * 2 / (n + 1); i++)
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

        /// <summary>
        /// Provides an Icon for the component.
        /// </summary>
        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                //You can add image files to your project resources and access them like this:
                // return Resources.IconForThisComponent;
                return null;
            }
        }

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("aa1e2686-32c8-43ab-a1f0-6b345dbd6ad4"); }
        }
    }
}