// Copyright (c) 2011 rubicon IT GmbH
using System;
using System.Collections.Generic;
using Emgu.CV;
using Emgu.CV.Structure;
using Emgu.CV.Util;
using Emgu.CV.XFeatures2D;
using System.Drawing;
using Emgu.CV.Flann;
using Emgu.CV.Features2D;
using Emgu.CV.CvEnum;
using System.Windows.Shapes;
using SW = System.Windows;

namespace UVAPositioning
{
    struct SIFTParametrs
    {
        public int nFeatures;
        public int nOctaveLayers;
        public double contrastThreshold;
        public double edgeThreshold;
        public double sigma;
    }
    class OriantatioOnMap
    {
        public Image<Rgb, byte> Map;
        public VectorOfKeyPoint VectorMapKeyPoint;
        public Mat MapDiscriptors = new Mat();
        public OriantatioOnMap(Image<Rgb, byte> Map, SIFTParametrs parametrs, double Compression = 4, double Radius = 20)
        {
            this.Map = Map;
            using (SIFT siftCPU = new SIFT(parametrs.nFeatures, parametrs.nOctaveLayers,
                 parametrs.contrastThreshold, parametrs.edgeThreshold, parametrs.sigma))
            {
                VectorMapKeyPoint = new VectorOfKeyPoint(siftCPU.Detect(Map));
                VectorMapKeyPoint = FilterKeyPoint(VectorMapKeyPoint, Map, Compression, Radius, parametrs);
                siftCPU.Compute(Map, VectorMapKeyPoint, MapDiscriptors);
            }
        }

        public Image<Rgb, byte> ShowMatches(Image<Rgb, byte> SubMap, int k, double uniquenessThreshold, 
            int gridx, int gridy, double persent, SIFTParametrs parametrs, out System.Windows.Point? location)
        {
            location = null;
            VectorOfKeyPoint VectorSubMapKeyPoint = null;
            Mat SubMapDiscriptors = null;
            VectorOfVectorOfDMatch matches = null;
            Mat mask = null;
            System.Drawing.Rectangle zone = new System.Drawing.Rectangle();
            Mat result = new Mat();
            Mat homography = null;
            try
            {

                FindMatches(SubMap, out VectorSubMapKeyPoint, out SubMapDiscriptors, out matches, out mask, out zone, out homography,
                    k, uniquenessThreshold, parametrs);
                Features2DToolbox.DrawMatches(SubMap, VectorSubMapKeyPoint, Map, VectorMapKeyPoint, matches,
                    result, new MCvScalar(0, 255, 0), new MCvScalar(0, 0, 255), mask, Features2DToolbox.KeypointDrawType.DrawRichKeypoints);
                PointF[] points = GetMapPoint(matches, mask);
                Point point = FoundCenter(points);
                if (MatchCorrect(Map.Mat, points, gridx, gridy, persent) && (!double.IsNaN(point.X) && !double.IsNaN(point.Y)))
                {
                    try
                    {
                        CvInvoke.Circle(result, point, 13, new MCvScalar(255, 0, 0), 10);
                        location = ImageTransform.SDtoSW(point);
                    }
                    catch { }
                }
                return new Image<Rgb, byte>(result.Bitmap);
            }
            catch
            {
                return Map;
            }
        }

        public void FindMatches(Image<Rgb, byte> SubMap, out VectorOfKeyPoint VectorSubMapKeyPoint,
            out Mat SubMapDiscriptors, out VectorOfVectorOfDMatch matches,
            out Mat mask, out System.Drawing.Rectangle zone, out Mat homography, int k, double uniquenessThreshold, SIFTParametrs parametrs)
        {
            VectorSubMapKeyPoint = new VectorOfKeyPoint();
            SubMapDiscriptors = new Mat();
            matches = new VectorOfVectorOfDMatch();
            zone = new System.Drawing.Rectangle();
            using (SIFT siftCPU = new SIFT(parametrs.nFeatures, parametrs.nOctaveLayers,
                parametrs.contrastThreshold, parametrs.edgeThreshold, parametrs.sigma))
            {
                siftCPU.DetectAndCompute(SubMap, null, VectorSubMapKeyPoint, SubMapDiscriptors, false);
            }
            matches = new VectorOfVectorOfDMatch();
            using (Emgu.CV.Flann.LinearIndexParams ip = new Emgu.CV.Flann.LinearIndexParams())
            using (Emgu.CV.Flann.SearchParams sp = new SearchParams())
            using (Emgu.CV.Features2D.DescriptorMatcher matcher = new FlannBasedMatcher(ip, sp))
            {
                matcher.Add(SubMapDiscriptors);
                matcher.KnnMatch(MapDiscriptors, matches, k, null);
            }

            mask = new Mat(matches.Size, 1, DepthType.Cv8U, 1);
            mask.SetTo(new MCvScalar(255));
            Features2DToolbox.VoteForUniqueness(matches, uniquenessThreshold, mask);

            homography = null;
            
            int nonZeroCount = CvInvoke.CountNonZero(mask);
            if (nonZeroCount >= 4)
            {
                nonZeroCount = Features2DToolbox.VoteForSizeAndOrientation(VectorSubMapKeyPoint, VectorMapKeyPoint,
                matches, mask, 1.5, 20);
                if (nonZeroCount >= 4)
                    homography = Features2DToolbox.GetHomographyMatrixFromMatchedFeatures(
                    VectorSubMapKeyPoint, VectorMapKeyPoint, matches, mask, 2);
            }

        }

        private VectorOfKeyPoint FilterKeyPoint(VectorOfKeyPoint InputVecor, Image<Rgb, byte> SourceImage, double Compression, double Diameter, SIFTParametrs parametrs)
        {
            VectorOfKeyPoint OutputVector = null;
            SourceImage = SourceImage.Resize(1.0 / Compression, Emgu.CV.CvEnum.Inter.Area);
            using (SIFT siftCPU = new SIFT(parametrs.nFeatures, parametrs.nOctaveLayers,
                parametrs.contrastThreshold, parametrs.edgeThreshold, parametrs.sigma))
            {
                VectorOfKeyPoint MainVecor = new VectorOfKeyPoint(siftCPU.Detect(SourceImage, null));
                OutputVector = new VectorOfKeyPoint(RemoveFakeKeyPoint(MainVecor, InputVecor, Compression, Diameter));
            }
            return OutputVector;
        }

        private MKeyPoint[] RemoveFakeKeyPoint(VectorOfKeyPoint MainVecor, VectorOfKeyPoint InputVecor, double Compression, double Radius)
        {
            List<MKeyPoint> InputListKeyPoint = new List<MKeyPoint>(InputVecor.ToArray());
            List<MKeyPoint> OutputVector = new List<MKeyPoint>();
            for (int i = 0; i < MainVecor.Size; i++)
            {
                for (int j = InputListKeyPoint.Count - 1; j >= 0; j--)
                {
                    PointF InputLocate = InputListKeyPoint[j].Point;
                    PointF MainLocate = MainVecor[i].Point;
                    if (Math.Pow(MainLocate.X * Compression - InputLocate.X, 2) +
                        Math.Pow(MainLocate.Y * Compression - InputLocate.Y, 2) <= Math.Pow(Radius, 2))
                    {
                        OutputVector.Add(InputListKeyPoint[j]);
                        InputListKeyPoint.RemoveAt(j);
                    }
                }
            }
            return OutputVector.ToArray();
        }

        private bool MatchCorrect(Mat image, PointF[] points, int gridx, int gridy, double persent)
        {
            double[,] f = new double[gridx, gridy];
            double sizex = (double)image.Width / gridx;
            double sizey = (double)image.Height/ gridy;
            for (int i = 0; i < points.Length; i++)
            {
                f[(int)(points[i].X / sizex), (int)(points[i].Y / sizey)] += 1.0 / points.Length;
            }
            double max = 0;
            for (int i = 0; i < gridx;i++)
            {
                for (int j = 0; j < gridy; j++)
                {
                    max = Math.Max(max, f[i, j]);
                }
            }
            return max >= persent;
        }
        
        private PointF[] GetMapPoint(VectorOfVectorOfDMatch matches, Mat mask)
        {
            List<PointF> points = new List<PointF>();
            for (int i = 0; i < matches.Size; i++)
            {
                var match = matches[i].ToArray();
                if (mask.GetData(i)[0] == 0)
                    continue;
                points.Add(VectorMapKeyPoint[match[0].QueryIdx].Point);
            }
            return points.ToArray();
        }

        private Point FoundCenter(PointF[] points)
        {
            PointF center = new PointF(0, 0);
            foreach (PointF point in points)
            {
                center.X += point.X;
                center.Y += point.Y;
            }
            center.X /= points.Length;
            center.Y /= points.Length;
            return new Point((int)Math.Round(center.X), (int)Math.Round(center.Y));
        }
    }
}
