using Emgu.CV;
using Emgu.CV.CvEnum;
using Emgu.CV.Features2D;
using Emgu.CV.Structure;
using System;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using SD = System.Drawing;

namespace UVAPositioning
{
    static class ImageTransform
    {


        static public Image<Rgb, byte> SetAircraft(Image<Rgb, byte> image, Aircraft aircraft)
        {
            Image<Rgb, byte> imageWithIcon = new Image<Rgb, byte>(image.Size);
            var rotateIcon = aircraft.aircraftIcon.Rotate(180 - aircraft.angle * 180 / Math.PI, new Rgb(Color.FromArgb(0, 0, 0, 0)));
            var resizeIcon = rotateIcon.Resize(0.08 * image.Height / aircraft.aircraftIcon.Height, Inter.Lanczos4);
            var roi = image.ROI;
            Bitmap imageWithIconBitmap = imageWithIcon.Bitmap;
            using (Graphics gr = Graphics.FromImage(imageWithIconBitmap))
            {
                SD.Point point = new Point()
                {
                    X = (int)(aircraft.CenterImage.X - resizeIcon.Width / 2.0),
                    Y = (int)(aircraft.CenterImage.Y - resizeIcon.Height / 2.0)
                };
                gr.DrawImage(resizeIcon.Bitmap, pjjoint);
            }
            imageWithIcon = new Image<Rgb, byte>(imageWithIconBitmap);
            return image + imageWithIcon;
        }
        static public Image<Rgb, byte> SetGrid(Image<Rgb, byte> image, System.Drawing.Point GridSize, Size size)
        {
            Image<Rgb, byte> result = image.Clone();
            for (int c = 1; c < GridSize.X; c++)
            {
                CvInvoke.Line(result,
                    new System.Drawing.Point(c * size.Width / GridSize.X, 0),
                    new System.Drawing.Point(c * size.Width / GridSize.X, size.Height),
                    new MCvScalar(125, 125, 125), 2);
            }
            for (int r = 1; r < GridSize.Y; r++)
            {
                CvInvoke.Line(result,
                    new System.Drawing.Point(0, r * size.Height / GridSize.Y),
                    new System.Drawing.Point(size.Width, r * size.Height / GridSize.Y),
                    new MCvScalar(125, 125, 125), 2);
            }
            return result;
        }

        static public Image<Rgb, byte> SetLine(Image<Rgb, byte> image, List<System.Windows.Point> WayAircraft, MCvScalar color)
        {
            Image<Rgb, byte> result = image.Clone();
            for (int i = 0; i < WayAircraft.Count - 1; i++)
            {
                System.Drawing.Point start = new System.Drawing.Point(Convert.ToInt32(WayAircraft[i].X),
                    Convert.ToInt32(WayAircraft[i].Y));
                System.Drawing.Point end = new System.Drawing.Point(Convert.ToInt32(WayAircraft[i + 1].X),
                    Convert.ToInt32(WayAircraft[i + 1].Y));
                CvInvoke.Line(result, start, end, color, 1, LineType.AntiAlias);
            }
            for (int i = 0; i < WayAircraft.Count; i++)
            {
                System.Drawing.Point element = new System.Drawing.Point(Convert.ToInt32(WayAircraft[i].X),
                    Convert.ToInt32(WayAircraft[i].Y));
                CvInvoke.Circle(result, element, 3, color, 2);
            }

            return result;
        }

        static public Image<Rgb, byte> SetKeyPont(Image<Rgb, byte> image, OriantatioOnMap oriantation, Bgr color)
        {
            Image<Rgb, byte> result = image.Clone();
            Features2DToolbox.DrawKeypoints(image, oriantation.VectorMapKeyPoint, result, color);
            return result;
        }
    }
}
