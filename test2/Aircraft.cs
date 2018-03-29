using Emgu.CV;
using Emgu.CV.Structure;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using SD = System.Drawing;

namespace UVAPositioning
{

    class Aircraft
    {
        List<System.Windows.Point> WayAircraft = new List<System.Windows.Point>();
        List<double> WaysLen = new List<double>();
        double WayLen;
        public double Locate = 0;
        public double Step = 0.05;
        public double Angle { get; private set; }
        public SD.Point CenterImage { get; private set; }

        public Image<Rgb, byte> AircraftIcon { get; private set; }

        public Aircraft(List<System.Windows.Point> WayAircraft, double Step, Image<Rgb, byte> aircraftIcon)
        {
            this.Step = Step;
            this.WayAircraft = WayAircraft;
            this.AircraftIcon = aircraftIcon;
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
            FindLocate();
        }
        public void FindLocate()
        {
            int index = 0;
            double LocateontheElement = 0;
            for (double l = 0; index < WaysLen.Count; l += WaysLen[index], index++)
                if (l + WaysLen[index] >= WayLen * Locate)
                {
                    LocateontheElement = (WayLen * Locate - l) / WaysLen[index];
                    break;
                }
            CenterImage =
                index != WaysLen.Count - 1 ?
                new SD.Point(
                (int)(LocateontheElement * (WayAircraft[index + 1].X - WayAircraft[index].X) + WayAircraft[index].X),
                (int)(LocateontheElement * (WayAircraft[index + 1].Y - WayAircraft[index].Y) + WayAircraft[index].Y)
                ) :
                new SD.Point(
                (int)(WayAircraft[index].X),
                (int)(WayAircraft[index].Y)
                );
            double y = WayAircraft[index + 1].Y - WayAircraft[index].Y;
            double x = WayAircraft[index + 1].X - WayAircraft[index].X;

            Angle = Math.Acos(y / Math.Sqrt(x * x + y * y));
            if (x < 0)
                Angle *= -1;
        }

        public Image<Rgb, byte> GetPhoto(Image<Rgb, byte> Map, SD.Point SizeOutImage)
        {

            Image<Rgb, byte> rotateImage = Map.Rotate(Angle * 180.0 / Math.PI, new Rgb(255, 255, 255), false);
            SD.Point centerRotateImage = new SD.Point()
            {
                X = (int)((CenterImage.X - Map.Width / 2.0) * Math.Cos(Angle) - 
                (CenterImage.Y - Map.Height / 2.0) * Math.Sin(Angle) + rotateImage.Width / 2.0),
                Y= (int)((CenterImage.X - Map.Width / 2.0) * Math.Sin(Angle) +
                (CenterImage.Y - Map.Height / 2.0) * Math.Cos(Angle) + rotateImage.Height/ 2.0)
            };

            var aircraftPhoto = new Image<Rgb, byte>(new Mat(rotateImage.Mat,
                new SD.Rectangle(
                    Math.Max(centerRotateImage.X - SizeOutImage.X / 2, 0),
                    Math.Max(centerRotateImage.Y - SizeOutImage.Y / 2, 0),
                    Math.Min(SizeOutImage.X, rotateImage.Mat.Width - centerRotateImage.X + SizeOutImage.X / 2 - 1),
                    Math.Min(SizeOutImage.Y, rotateImage.Mat.Height - centerRotateImage.Y + SizeOutImage.Y / 2 - 1)
                    )).Bitmap
                    );
            Locate += Step;
            return aircraftPhoto.Rotate(180,new Rgb(255,255,255));
        }

        private SD.Point RotatePoint(SD.Point point, SD.Point center, double angle)
        {
            SD.Point result = new SD.Point()
            {
                X = (int)((point.X - center.X) * Math.Cos(angle) - (point.X - center.X) * Math.Sin(angle)) + center.X,
                Y = (int)((point.X - center.X) * Math.Sin(angle) + (point.X - center.X) * Math.Cos(angle)) + center.Y
            };
            return result;
        }
    }
}
