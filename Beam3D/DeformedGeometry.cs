using System;
using System.Collections.Generic;

using Grasshopper.Kernel;
using Rhino.Geometry;

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


        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddNumberParameter("Deformation", "Def", "Deformations from 3DBeamCalc", GH_ParamAccess.list);
            pManager.AddLineParameter("Geometry", "G", "Input Geometry (Line format)", GH_ParamAccess.list);
            pManager.AddNumberParameter("Scale", "S", "The Scale Factor for Deformation", GH_ParamAccess.item, 1000);
        }

        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddLineParameter("Deformed Geometry", "Def.G.", "Deformed Geometry as List of Lines", GH_ParamAccess.list);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            #region Fetch input
            //Expected inputs and outputs
            List<double> def = new List<double>();
            List<Line> geometry = new List<Line>();
            double scale = 1000;
            List<Line> defGeometry = new List<Line>();
            List<Point3d> defPoints = new List<Point3d>();

            //Set expected inputs from Indata
            if (!DA.GetDataList(0, def)) return;
            if (!DA.GetDataList(1, geometry)) return;
            if (!DA.GetData(2, ref scale)) return;
            #endregion

            #region Scale deformations
            //List all nodes (every node only once), numbering them according to list index
            List<Point3d> points = CreatePointList(geometry);

            int index = 0;
            //loops through all points and scales x-, y- and z-dir
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
            #endregion

            #region Create geometry
            //creates deformed geometry based on initial geometry placement
            foreach (Line line in geometry)
            {
                //fetches index of original start and endpoint
                int i1 = points.IndexOf(line.From);
                int i2 = points.IndexOf(line.To);

                //creates new line based on scaled deformation of said points
                defGeometry.Add(new Line(defPoints[i1], defPoints[i2]));
            }
            #endregion


            //Set output data
            DA.SetDataList(0, defGeometry);
        }   //End of main program


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
    }
}