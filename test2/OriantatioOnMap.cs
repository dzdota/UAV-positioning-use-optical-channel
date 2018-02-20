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
        public OriantatioOnMap(Image<Rgb, byte> Map, SIFTParametrs parametrs, double Compression = 4, double Diameter = 20)
        {
            
            this.Map = Map;
            using (SIFT siftCPU = new SIFT(parametrs.nFeatures, parametrs.nOctaveLayers,
                 parametrs.contrastThreshold, parametrs.edgeThreshold, parametrs.sigma))
            {
                VectorMapKeyPoint = new VectorOfKeyPoint(siftCPU.Detect(Map));
                VectorMapKeyPoint = FilterKeyPoint(VectorMapKeyPoint, Map, Compression, Diameter, parametrs);
                siftCPU.Compute(Map, VectorMapKeyPoint, MapDiscriptors);
            }
        }

        public Image<Rgb, byte> ShowMatches(Image<Rgb, byte> SubMap, int k, double uniquenessThreshold, SIFTParametrs parametrs)
        {
            VectorOfKeyPoint VectorSubMapKeyPoint = null;
            Mat SubMapDiscriptors = null;
            VectorOfVectorOfDMatch matches = null;
            Mat mask = null;
            Rectangle zone = new Rectangle();
            Mat result = new Mat();
            FindLocateMap(SubMap, ref VectorSubMapKeyPoint, ref SubMapDiscriptors, ref matches, ref mask, ref zone, k, uniquenessThreshold, parametrs);
            Features2DToolbox.DrawMatches(SubMap, VectorSubMapKeyPoint, Map, VectorMapKeyPoint, matches, 
                result, new MCvScalar(0, 255, 0), new MCvScalar(0, 0, 255), mask, Features2DToolbox.KeypointDrawType.DrawRichKeypoints);
            return new Image<Rgb, byte>(result.Bitmap);
        }

        public void FindLocateMap(Image<Rgb, byte> SubMap, ref VectorOfKeyPoint VectorSubMapKeyPoint, 
            ref Mat SubMapDiscriptors, ref VectorOfVectorOfDMatch matches, ref Mat mask, ref Rectangle zone, int k, double uniquenessThreshold, SIFTParametrs parametrs)
        {
            VectorSubMapKeyPoint = new VectorOfKeyPoint();
            SubMapDiscriptors = new Mat();
            matches = new VectorOfVectorOfDMatch();
            zone = new Rectangle();
            using (SIFT siftCPU = new SIFT(parametrs.nFeatures, parametrs.nOctaveLayers,
                parametrs.contrastThreshold, parametrs.edgeThreshold, parametrs.sigma))
                siftCPU.DetectAndCompute(SubMap, null, VectorSubMapKeyPoint, SubMapDiscriptors, false);
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

            Mat homography = null;

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
        private MKeyPoint[] RemoveFakeKeyPoint(VectorOfKeyPoint MainVecor, VectorOfKeyPoint InputVecor, double Compression, double Diameter)
        {
            List<MKeyPoint> InputListKeyPoint = new List<MKeyPoint>(InputVecor.ToArray());
            List<MKeyPoint> OutputVector = new List<MKeyPoint>();
            for (int i = 0; i < MainVecor.Size; i++)
            {
                for (int j = InputListKeyPoint.Count - 1; j >= 0; j--)
                {
                    PointF InputLocate = InputListKeyPoint[j].Point;
                    PointF MainLocate = MainVecor[i].Point;
                    if (Math.Pow(MainLocate.X * Compression - InputLocate.X, 2) + Math.Pow(MainLocate.Y * Compression - InputLocate.Y, 2) <= Math.Pow(Diameter, 2))
                    {
                        OutputVector.Add(InputListKeyPoint[j]);
                        InputListKeyPoint.RemoveAt(j);
                    }
                }
            }
            return OutputVector.ToArray();
        }
    }
}
