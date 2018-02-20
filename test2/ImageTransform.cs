using Emgu.CV;
using Emgu.CV.CvEnum;
using Emgu.CV.Features2D;
using Emgu.CV.Structure;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace test2
{
    static class ImageTransform
    {


        static public Image<Rgb, byte> SetGrid(Image<Rgb, byte> image, System.Drawing.Point GridSize)
        {
            Image<Rgb, byte> result = image.Clone();
            for (int c = 1; c < GridSize.X; c++)
            {
                CvInvoke.Line(result,
                    new System.Drawing.Point(c * result.Width / GridSize.X, 0),
                    new System.Drawing.Point(c * result.Width / GridSize.X, result.Height),
                    new MCvScalar(125, 125, 125), 2);
            }
            for (int r = 1; r < GridSize.Y; r++)
            {
                CvInvoke.Line(result,
                    new System.Drawing.Point(0, r * result.Height / GridSize.Y),
                    new System.Drawing.Point(result.Width, r * result.Height / GridSize.Y),
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
