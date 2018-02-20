using Emgu.CV;
using Emgu.CV.Structure;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace test2
{

    class Aircraft
    {
        List<System.Windows.Point> WayAircraft = new List<System.Windows.Point>();
        List<double> WaysLen = new List<double>();
        double WayLen;
        public double Locate = 0;
        public double Step = 0.1;
        public Aircraft(List<System.Windows.Point> WayAircraft, double Step)
        {
            this.Step = 0.1;
            this.WayAircraft = WayAircraft;
            if (Step > 0 && Step < 1)
                this.Step = Step;
            WayLen = 0;
            WaysLen.Clear();
            for (int i = 0; i < WayAircraft.Count - 1; i++)
            {
                WaysLen.Add(Math.Sqrt(
                    Math.Pow(WayAircraft[i + 1].X - WayAircraft[i].X, 2) +
                    Math.Pow(WayAircraft[i + 1].Y - WayAircraft[i].Y, 2)
                    ));
                WayLen += WaysLen[WaysLen.Count - 1]; 
            }
        }

        public Image<Rgb, byte> NextLocate(Image<Rgb, byte> Map, System.Drawing.Point SizeOutImage)
        {
            //Mat res = new Mat(Map.Mat, )
            //Image<Rgb, byte> result = new Image<Rgb, byte>(;
            int index = 0;
            double LocateontheElement = 0;
            for (double l = 0; index < WaysLen.Count; l += WaysLen[index], index++)  
                if (l + WaysLen[index] >= WayLen * Locate)
                {
                    LocateontheElement = (WayLen * Locate - l) / WaysLen[index];
                    break;
                }
            System.Drawing.Point CenterImage = 
                LocateontheElement != 0 ?
                new System.Drawing.Point(
                (int)(LocateontheElement * (WayAircraft[index + 1].X - WayAircraft[index].X) + WayAircraft[index].X),
                (int)(LocateontheElement * (WayAircraft[index + 1].Y - WayAircraft[index].Y) + WayAircraft[index].Y)
                ) :
                new System.Drawing.Point(
                (int)(WayAircraft[index].X),
                (int)(WayAircraft[index].Y)
                );/*
            double angle = Math.Atan((WayAircraft[index + 1].Y - WayAircraft[index].Y) / (WayAircraft[index + 1].X - WayAircraft[index].X));
            System.Drawing.Point UpperLeft = RotatePoint(
                new System.Drawing.Point(CenterImage.X - SizeOutImage.X / 2, CenterImage.Y - SizeOutImage.Y / 2), CenterImage, angle);
            System.Drawing.Point UpperRight = RotatePoint(
                new System.Drawing.Point(CenterImage.X + SizeOutImage.X / 2, CenterImage.Y - SizeOutImage.Y / 2), CenterImage, angle);
            System.Drawing.Point BottomLeft = RotatePoint(
                new System.Drawing.Point(CenterImage.X - SizeOutImage.X / 2, CenterImage.Y + SizeOutImage.Y / 2), CenterImage, angle);
            System.Drawing.Point BottomRight = RotatePoint(
                new System.Drawing.Point(CenterImage.X + SizeOutImage.X / 2, CenterImage.Y + SizeOutImage.Y / 2), CenterImage, angle);*/

            Locate += Step;

            return new Image<Rgb, byte>( new Mat(Map.Mat,
                new System.Drawing.Rectangle(
                    Math.Max(CenterImage.X - SizeOutImage.X / 2, 0),
                    Math.Max(CenterImage.Y - SizeOutImage.Y / 2, 0),
                    Math.Min(SizeOutImage.X, Map.Mat.Width -  CenterImage.X + SizeOutImage.X / 2 - 1),
                    Math.Min(SizeOutImage.Y, Map.Mat.Height - CenterImage.Y + SizeOutImage.Y / 2 - 1)
                    )).Bitmap
                    );
        }

        private System.Drawing.Point RotatePoint(System.Drawing.Point point, System.Drawing.Point center, double angle)
        {
            System.Drawing.Point result = new System.Drawing.Point()
            {
                X = (int)((point.X - center.X) * Math.Cos(angle) - (point.X - center.X) * Math.Sin(angle)) + center.X,
                Y = (int)((point.X - center.X) * Math.Sin(angle) + (point.X - center.X) * Math.Cos(angle)) + center.Y
            };
            return result;
        }
    }
}
