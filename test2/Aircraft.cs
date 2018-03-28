using Emgu.CV;
using Emgu.CV.Structure;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace UVAPositioning
{

    class Aircraft
    {
        List<System.Windows.Point> WayAircraft = new List<System.Windows.Point>();
        List<double> WaysLen = new List<double>();
        double WayLen;
        public double Locate = 0;
        public double Step = 0.05;
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
        public System.Drawing.Point FindLocate(out double angle)
        {
            int index = 0;
            double LocateontheElement = 0;
            for (double l = 0; index < WaysLen.Count; l += WaysLen[index], index++)
                if (l + WaysLen[index] >= WayLen * Locate)
                {
                    LocateontheElement = (WayLen * Locate - l) / WaysLen[index];
                    break;
                }
            System.Drawing.Point CenterImage =
                index != WaysLen.Count - 1 ?
                new System.Drawing.Point(
                (int)(LocateontheElement * (WayAircraft[index + 1].X - WayAircraft[index].X) + WayAircraft[index].X),
                (int)(LocateontheElement * (WayAircraft[index + 1].Y - WayAircraft[index].Y) + WayAircraft[index].Y)
                ) :
                new System.Drawing.Point(
                (int)(WayAircraft[index].X),
                (int)(WayAircraft[index].Y)
                );
            double y = WayAircraft[index + 1].Y - WayAircraft[index].Y;
            double x = WayAircraft[index + 1].X - WayAircraft[index].X;

            angle = Math.Acos(y / Math.Sqrt(x * x + y * y));
            return CenterImage;
        }

        public Image<Rgb, byte> GetPhoto(Image<Rgb, byte> Map, System.Drawing.Point SizeOutImage)
        {
            double angle;
            System.Drawing.Point CenterImage = FindLocate(out angle);

            Image<Rgb, byte> rotateImage = Map.Rotate(angle * 180.0 / Math.PI, new Rgb(255, 255, 255), false);
            System.Drawing.Point centerRotateImage = new System.Drawing.Point()
            {
                X = (int)((CenterImage.X - Map.Width / 2.0) * Math.Cos(angle) - 
                (CenterImage.Y - Map.Height / 2.0) * Math.Sin(angle) + rotateImage.Width / 2.0),
                Y= (int)((CenterImage.X - Map.Width / 2.0) * Math.Sin(angle) +
                (CenterImage.Y - Map.Height / 2.0) * Math.Cos(angle) + rotateImage.Height/ 2.0)
            };

            var aircraftPhoto = new Image<Rgb, byte>(new Mat(rotateImage.Mat,
                new System.Drawing.Rectangle(
                    Math.Max(centerRotateImage.X - SizeOutImage.X / 2, 0),
                    Math.Max(centerRotateImage.Y - SizeOutImage.Y / 2, 0),
                    Math.Min(SizeOutImage.X, rotateImage.Mat.Width - centerRotateImage.X + SizeOutImage.X / 2 - 1),
                    Math.Min(SizeOutImage.Y, rotateImage.Mat.Height - centerRotateImage.Y + SizeOutImage.Y / 2 - 1)
                    )).Bitmap
                    );
            Locate += Step;
            return aircraftPhoto.Rotate(180,new Rgb(255,255,255));
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
