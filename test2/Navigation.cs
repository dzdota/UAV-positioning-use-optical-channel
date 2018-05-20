using Emgu.CV;
using Emgu.CV.Structure;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using SD = System.Drawing;
using SW = System.Windows;

namespace UVAPositioning
{
    public class Icon
    {
        public Image<Rgb, byte> IconImage { get; private set; }
        public SW.Point Coordinate { get; private set; }
        public Icon(Image<Rgb, byte> iconImage, SW.Point coordinate)
        {
            this.IconImage = iconImage;
            this.Coordinate = coordinate;
        }
    }
    class Navigation
    {
        Icon baseIcon;
        Icon lostSignalIcon;
        public Point? location;
        public Point? Location {
            get
            {
                return location;
            }
            set
            {
                location = value;
                sure = value == null ? 0 : 1;
            }
        }
        Point lost;
        Point cameraSize;
        OriantatioOnMap OoM;
        double angle = 0;
        double k;
        double sure = 0;

        public Navigation(OriantatioOnMap OoM, Icon baseIcon, Icon lostSignalIcon, Point cameraSize, Point? location = null)
        {
            this.OoM = OoM;
            this.baseIcon = baseIcon;
            this.lostSignalIcon = lostSignalIcon;
            this.lost = lostSignalIcon.Coordinate;
            this.cameraSize = cameraSize;
            this.k = cameraSize.X * 0.95 / (2 * Math.PI);
            this.Location = location;
        }

        public SW.Point NextLocation()
        {
            double step = cameraSize.Y * 0.95;
            if (Location == null)
            {
                double L = LenSpiral(angle);
                angle = AngleSpiralIterate(L + step);
                return new SW.Point(lostSignalIcon.Coordinate.X +  Math.Cos(angle) * angle * k,
                    lostSignalIcon.Coordinate.Y + Math.Sin(angle) * angle * k);
            }
            else
            {
                double Len = LenLine(location.Value, baseIcon.Coordinate);
                location = new Point(
                    location.Value.X + (step / Len) * (baseIcon.Coordinate.X - location.Value.X),
                    location.Value.Y + (step / Len) * (baseIcon.Coordinate.Y - location.Value.Y));
                sure -= 0.05;
                return location.Value;
            }
        }
        private double LenLine(Point p1, Point p2)
        {
            return Math.Sqrt(Math.Pow(p1.X - p2.X, 2) + Math.Pow(p1.Y - p2.Y, 2));
        }

        private double AngleSpiralIterate(double Len)
        {
            double e = 0.00001;
            double angle = 0.002, anglePrev = -1;
            while (Math.Abs(angle - anglePrev) >e)
            {
                anglePrev = angle;
                angle = angle -  (LenSpiral(angle) - Len) / (Math.Sqrt(1 + angle * angle) * k);
            }
            return angle;
        }
        private double  LenSpiral( double angle)
        {
            return (k / 2) * 
                (angle * Math.Sqrt(1 + angle * angle) 
                + Math.Log(
                    angle + Math.Sqrt(1 + angle * angle)
                    ));
        }

        public Image<Rgb, byte> GetPhoto(Image<Rgb, byte> Map, SW.Point locationPhoto, SW.Point SizeOutImage)
        {
            double angleRotate = Math.Atan((Math.Sqrt(1 + angle * angle) * cameraSize.X * 0.95));
            Image<Rgb, byte> rotateImage = Map.Rotate(angleRotate * 180.0 / Math.PI, new Rgb(255, 255, 255), false);
            SW.Point centerRotateImage = new SW.Point()
            {
                X = (int)((locationPhoto.X - Map.Width / 2.0) * Math.Cos(angleRotate) -
                (locationPhoto.Y - Map.Height / 2.0) * Math.Sin(angleRotate) + rotateImage.Width / 2.0),
                Y = (int)((locationPhoto.X - Map.Width / 2.0) * Math.Sin(angleRotate) +
                (locationPhoto.Y - Map.Height / 2.0) * Math.Cos(angleRotate) + rotateImage.Height / 2.0)
            };

            var aircraftPhoto = new Image<Rgb, byte>(new Mat(rotateImage.Mat,
                new SD.Rectangle(
                    (int)Math.Max(centerRotateImage.X - SizeOutImage.X / 2, 0),
                    (int)Math.Max(centerRotateImage.Y - SizeOutImage.Y / 2, 0),
                    (int)Math.Min(SizeOutImage.X, rotateImage.Mat.Width - centerRotateImage.X + SizeOutImage.X / 2 - 1),
                    (int)Math.Min(SizeOutImage.Y, rotateImage.Mat.Height - centerRotateImage.Y + SizeOutImage.Y / 2 - 1)
                    )).Bitmap
                    );
            return aircraftPhoto.Rotate(180, new Rgb(255, 255, 255));
        }
    }
}
