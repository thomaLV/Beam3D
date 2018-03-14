using System;
using System.Collections.Generic;

using Grasshopper.Kernel;
using Rhino.Geometry;

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

        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddPointParameter("Points", "P", "Points to apply moment", GH_ParamAccess.list);
            pManager.AddNumberParameter("Moment", "M", "Moment Magnitude [kNm]", GH_ParamAccess.list);
            pManager.AddBooleanParameter("Mx", "axz", "Add Moment around X-axis?", GH_ParamAccess.list);
            pManager.AddBooleanParameter("My", "axy", "Add Moment around Y-axis?", GH_ParamAccess.list);
            pManager.AddBooleanParameter("Mz", "axy", "Add Moment around Z-axis? ", GH_ParamAccess.list);
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
            List<bool> mx = new List<bool>();
            List<bool> my = new List<bool>();
            List<bool> mz = new List<bool>();
            List<string> pointInStringFormat = new List<string>();      //preallocate final string output

            //Set expected inputs from Indata
            if (!DA.GetDataList(0, pointList)) return;
            if (!DA.GetDataList(1, momentList)) return;
            DA.GetDataList(2, mx);
            DA.GetDataList(3, my);
            DA.GetDataList(4, mz);
            #endregion

            #region Prepare variables
            //initialize temporary stringline and moment magnitudes
            string mxs = "0";
            string mys = "0";
            string mzs = "0";
            string vectorString;

            //populate missing entries in momentlist with copies of last entry
            if (momentList.Count < pointList.Count)
            {
                for (int i = pointList.Count - (pointList.Count - momentList.Count); i < pointList.Count; i++)
                {
                    momentList.Add(momentList[momentList.Count - 1]);
                }
            }
            #endregion

            #region Format output if mx, my or mz has special considerations
            if (mx.Count == 1 && my.Count != 1 && mz.Count != 1)
            {
                for (int i = 0; i < pointList.Count; i++)
                {
                    string magnitude = momentList[i].ToString();
                    if (mx[0]) mxs = magnitude;
                    if (my[i])
                    {
                        mys = magnitude;
                    }
                    else
                    {
                        mys = "0";
                    }
                    if (mz[i])
                    {
                        mzs = magnitude;
                    }
                    else
                    {
                        mzs = "0";
                    }

                    vectorString = mxs + "," + mys + "," + mzs;
                    pointInStringFormat.Add(pointList[i].X + "," + pointList[i].Y + "," + pointList[i].Z + ":" + vectorString);
                }
            }
            else if (mx.Count != 1 && my.Count == 1 && mz.Count != 1)
            {
                for (int i = 0; i < pointList.Count; i++)
                {
                    string magnitude = momentList[i].ToString();
                    if (mx[i])
                    {
                        mxs = magnitude;
                    }
                    else
                    {
                        mxs = "0";
                    }
                    if (my[0]) mys = magnitude;
                    if (mz[i])
                    {
                        mzs = magnitude;
                    }
                    else
                    {
                        mzs = "0";
                    }

                    vectorString = mxs + "," + mys + "," + mzs;
                    pointInStringFormat.Add(pointList[i].X + "," + pointList[i].Y + "," + pointList[i].Z + ":" + vectorString);
                }
            }
            else if (mx.Count != 1 && my.Count != 1 && mz.Count == 1)
            {
                for (int i = 0; i < pointList.Count; i++)
                {
                    string magnitude = momentList[i].ToString();
                    if (mx[i])
                    {
                        mxs = magnitude;
                    }
                    else
                    {
                        mxs = "0";
                    }
                    if (my[i])
                    {
                        mys = magnitude;
                    }
                    else
                    {
                        mys = "0";
                    }
                    if (mz[0]) mzs = magnitude;

                    vectorString = mxs + "," + mys + "," + mzs;
                    pointInStringFormat.Add(pointList[i].X + "," + pointList[i].Y + "," + pointList[i].Z + ":" + vectorString);
                }
            }
            else if (mx.Count == 1 && my.Count == 1 && mz.Count != 1)
            {
                for (int i = 0; i < pointList.Count; i++)
                {
                    string magnitude = momentList[i].ToString();
                    if (mx[0]) mxs = magnitude;
                    if (my[i]) mys = magnitude;
                    if (mz[i])
                    {
                        mzs = magnitude;
                    }
                    else
                    {
                        mzs = "0";
                    }

                    vectorString = mxs + "," + mys + "," + mzs;
                    pointInStringFormat.Add(pointList[i].X + "," + pointList[i].Y + "," + pointList[i].Z + ":" + vectorString);
                }
            }
            else if (mx.Count == 1 && my.Count != 1 && mz.Count == 1)
            {
                for (int i = 0; i < pointList.Count; i++)
                {
                    string magnitude = momentList[i].ToString();
                    if (mx[0]) mxs = magnitude;
                    if (my[i])
                    {
                        mys = magnitude;
                    }
                    else
                    {
                        mys = "0";
                    }
                    if (mz[0]) mzs = magnitude;


                    vectorString = mxs + "," + mys + "," + mzs;
                    pointInStringFormat.Add(pointList[i].X + "," + pointList[i].Y + "," + pointList[i].Z + ":" + vectorString);
                }
            }
            else if (mx.Count != 1 && my.Count == 1 && mz.Count == 1)
            {
                for (int i = 0; i < pointList.Count; i++)
                {
                    string magnitude = momentList[i].ToString();
                    if (mx[i])
                    {
                        mxs = magnitude;
                    }
                    else
                    {
                        mxs = "0";
                    }
                    if (my[0]) mys = magnitude;
                    if (mz[0]) mzs = magnitude;

                    vectorString = mxs + "," + mys + "," + mzs;
                    pointInStringFormat.Add(pointList[i].X + "," + pointList[i].Y + "," + pointList[i].Z + ":" + vectorString);
                }
            }
            else if (mx.Count == 1 && my.Count == 1 && mz.Count == 1)
            {
                for (int i = 0; i < pointList.Count; i++)
                {
                    string magnitude = momentList[i].ToString();
                    if (mx[0]) mxs = magnitude;
                    if (my[0]) mys = magnitude;
                    if (mz[0]) mzs = magnitude;

                    vectorString = mxs + "," + mys + "," + mzs;
                    pointInStringFormat.Add(pointList[i].X + "," + pointList[i].Y + "," + pointList[i].Z + ":" + vectorString);
                }
            }
            #endregion

            #region Format output otherwise
            else
            {
                for (int i = 0; i < pointList.Count; i++)
                {
                    string magnitude = momentList[i].ToString();
                    if (mx[i])
                    {
                        mxs = magnitude;
                    }
                    else
                    {
                        mxs = "0";
                    }
                    if (my[i])
                    {
                        mys = magnitude;
                    }
                    else
                    {
                        mys = "0";
                    }
                    if (mz[i])
                    {
                        mzs = magnitude;
                    }
                    else
                    {
                        mzs = "0";
                    }

                    vectorString = mxs + "," + mys + "," + mzs;
                    pointInStringFormat.Add(pointList[i].X + "," + pointList[i].Y + "," + pointList[i].Z + ":" + vectorString);
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
    }
}