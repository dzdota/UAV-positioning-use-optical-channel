
using System;
using System.Drawing;
using System.Windows;
using System.Windows.Input;
using System.Windows.Media.Imaging;
using System.IO;

using Emgu.CV;
using Emgu.CV.Structure;
using Emgu.CV.XFeatures2D;
using Emgu.CV.Util;
using Emgu.CV.Flann;
using Emgu.CV.Features2D;
using Emgu.CV.CvEnum;
using Emgu.CV.Cuda;
using Emgu.CV.ML;

using System.Collections.Generic;
using System.Windows.Forms;


using OpenCL.Net;

namespace test2
{
    /// <summary>
    /// Interaction logic for MainWindow.xaml
    /// </summary>
    public partial class MainWindow : Window
    {
        List<System.Windows.Point> WayAircraft = new List<System.Windows.Point>();
        OriantatioOnMap oriantation;
        Aircraft aircraft = null;
        Image<Rgb, byte> imagefromAircraft = null;
        SIFTParametrs parametrs;
        System.Drawing.Point GridSize;
        double Compression;
        double Diameter;
        public MainWindow()
        {
            InitializeComponent();
            parametrs = new SIFTParametrs()
            {
                nFeatures = 0,
                nOctaveLayers = 5,
                contrastThreshold = 0.1,
                edgeThreshold = 10,
                sigma = 1.6
            };
            
        }

        public void SetImageAircraft(Image<Rgb, byte> image)
        {
            /*double Transform = Math.Max(oriantation.Map.Width / imageBox1.Width,
                oriantation.Map.Height / imageBox1.Height);
            image = image.Resize(Transform, Inter.Area);*/
            image = ImageTransform.SetLine(image, WayAircraft, new MCvScalar(0, 0, 255));
            var imgbrush = new BitmapImage();
            imgbrush.BeginInit();
            var ms = new MemoryStream();
            image.Bitmap.Save(ms, System.Drawing.Imaging.ImageFormat.Bmp);
            imgbrush.StreamSource = ms;
            imgbrush.CreateOptions = BitmapCreateOptions.PreservePixelFormat;
            imgbrush.EndInit();
            imageBoxfromAircraft.Source = imgbrush;
        }
        public void SetImage(Image<Rgb, byte> image)
        {
            /*double Transform = Math.Max(oriantation.Map.Width / imageBox1.Width,
                oriantation.Map.Height / imageBox1.Height);
            image = image.Resize(Transform, Inter.Area);*/
            if (ShowGridCheckBox.IsChecked == true)
                image = ImageTransform.SetGrid(image, GridSize);
            if (ShowKeyPointCheckBox.IsChecked == true && oriantation != null)
                image = ImageTransform.SetKeyPont(image, oriantation, new Bgr(255, 0, 0));
            var imgbrush = new BitmapImage();
            imgbrush.BeginInit();
            var ms = new MemoryStream();
            image.Bitmap.Save(ms, System.Drawing.Imaging.ImageFormat.Bmp);
            imgbrush.StreamSource = ms;
            imgbrush.CreateOptions = BitmapCreateOptions.PreservePixelFormat;
            imgbrush.EndInit();
            imageBox1.Source = imgbrush;
            image.Save("outputimage.bmp");
        }

        [Obsolete("")]
        public Bitmap drawSift(Image<Rgb, byte> modelimage, Image<Rgb, byte> observedimage)
        {
            int k = 2;
            double uniquenessThreshold = 0.80;
            

            VectorOfKeyPoint modelKeyPoints = new VectorOfKeyPoint(),
                observedKeyPoints = new VectorOfKeyPoint();

            Mat modeldiscriptors = new Mat();
            Mat observeddiscriptors = new Mat();
            //observedKeyPoints = observedKeyPoints.Resize(1.0 / Compression, Inter.Area);
            using (SIFT siftCPU = new SIFT(0, 5, 0.04, 10.0, 1.6))
            {
                siftCPU.DetectAndCompute(modelimage, null, modelKeyPoints, modeldiscriptors, false);
                observedKeyPoints = new VectorOfKeyPoint(siftCPU.Detect(observedimage));
                siftCPU.Compute(observedimage, observedKeyPoints, observeddiscriptors);
            }
            VectorOfVectorOfDMatch matches = new VectorOfVectorOfDMatch();

            using (Emgu.CV.Flann.LinearIndexParams ip = new Emgu.CV.Flann.LinearIndexParams())
            using (Emgu.CV.Flann.SearchParams sp = new SearchParams())
            using (Emgu.CV.Features2D.DescriptorMatcher matcher = new FlannBasedMatcher(ip, sp))
            {
                matcher.Add(modeldiscriptors);
                matcher.KnnMatch(observeddiscriptors, matches, k, null);
            }

            Mat mask = new Mat(matches.Size, 1, DepthType.Cv8U, 1);
            mask.SetTo(new MCvScalar(255));
            Features2DToolbox.VoteForUniqueness(matches, uniquenessThreshold, mask);

            Mat homography = null;

            int nonZeroCount = CvInvoke.CountNonZero(mask);
            if (nonZeroCount >= 4)
            {
                nonZeroCount = Features2DToolbox.VoteForSizeAndOrientation(modelKeyPoints, observedKeyPoints,
                matches, mask, 1.5, 20);
                if (nonZeroCount >= 4)
                    homography = Features2DToolbox.GetHomographyMatrixFromMatchedFeatures(modelKeyPoints,
                    observedKeyPoints, matches, mask, 2);
            }

            observedimage = new Image<Rgb, byte>(DrawZone(observedimage.Mat, observedKeyPoints, matches, mask).Bitmap);
            //modelKeyPoints.FilterByPixelsMask(new Image<Gray, byte>(mask.Bitmap));
            //observedKeyPoints.FilterByPixelsMask(new Image<Gray, byte>(mask.Bitmap));

            Mat result = new Mat();
            //Draw the matched keypoints
            Features2DToolbox.DrawMatches(modelimage, modelKeyPoints, observedimage, observedKeyPoints,
            matches, result, new MCvScalar(0, 255, 0), new MCvScalar(255,0, 0), mask);
            if (homography != null)
            {
                //draw a rectangle along the projected model
                System.Drawing.Rectangle rect = new System.Drawing.Rectangle(System.Drawing.Point.Empty, modelimage.Size);
                PointF[] pts = new PointF[]
                {
                    new PointF(rect.Left, rect.Bottom),
                    new PointF(rect.Right, rect.Bottom),
                    new PointF(rect.Right, rect.Top),
                    new PointF(rect.Left, rect.Top)
                };
                pts = CvInvoke.PerspectiveTransform(pts, homography);
                
#if NETFX_CORE
                Point[] points = Extensions.ConvertAll<PointF, Point>(pts, Point.Round);
#else
                System.Drawing.Point[] points = Array.ConvertAll<PointF, System.Drawing.Point>(pts, System.Drawing.Point.Round);
#endif
                using (VectorOfPoint vp = new VectorOfPoint(points))
                {
                    CvInvoke.Polylines(result, vp, true, new MCvScalar(0, 0, 255), 2);
                }
            }

            return result.Bitmap;

        }

        public Mat DrawZone(Mat observedimage, VectorOfKeyPoint observedKeyPoints, VectorOfVectorOfDMatch matches, Mat mask, int n = 4)
        {
            Mat result = observedimage.Clone();
            int ypart = 10, xpart = 16;
            int[,] matrix = new int[ypart, xpart];
            for (int i = 0; i < matches.Size; i++)
            {
                var a = matches[i].ToArray();
                if (mask.GetData(i)[0] == 0)
                    continue;
                System.Drawing.Point p =
                    new System.Drawing.Point((int)observedKeyPoints[a[0].QueryIdx].Point.X,
                        (int)observedKeyPoints[a[0].QueryIdx].Point.Y);
                matrix[(int)p.Y * ypart / observedimage.Bitmap.Height, 
                    (int)p.X * xpart / observedimage.Bitmap.Width]++;
            }
            for (int x = 0; x < xpart; x++)
                for (int y = 0; y < ypart; y++)
                {
                    if (matrix[y, x] >= n)
                    {
                        System.Drawing.Rectangle rec = new System.Drawing.Rectangle(
                            (int)x * observedimage.Bitmap.Width / xpart,
                            (int)y * observedimage.Bitmap.Height / ypart,
                            observedimage.Bitmap.Width / xpart,
                            observedimage.Bitmap.Height / ypart);
                        CvInvoke.Rectangle(result, rec, new MCvScalar(0 ,0 , 255), 2);
                    }
                }
            return result;
        }

        private void imageBox1_MouseRightButtonDown(object sender, MouseButtonEventArgs e)
        {
            System.Windows.Point position = e.GetPosition(imageBoxfromAircraft);
            double Transform = Math.Max(imagefromAircraft.Width / imageBoxfromAircraft.Width,
                imagefromAircraft.Height / imageBoxfromAircraft.Height);
            position = new System.Windows.Point(position.X * Transform,
                position.Y * Transform);
            WayAircraft.Add(position);
            SetImageAircraft(imagefromAircraft);
        }

        private void ReadAircraftImage_Click(object sender, RoutedEventArgs e)
        {
            WayAircraft.Clear();
            OpenFileDialog openFileDialog = new OpenFileDialog();
            if (openFileDialog.ShowDialog() == System.Windows.Forms.DialogResult.OK)
            {
                try
                {
                    Mat SourceMat = new Mat(openFileDialog.FileName, ImreadModes.AnyColor);
                    imagefromAircraft = new Image<Rgb, byte>(SourceMat.Bitmap);
                    SetImageAircraft(imagefromAircraft);
                }
                catch { }
            }
        }

        private void OpenFile_Click(object sender, RoutedEventArgs e)
        {
            OpenFileDialog openFileDialog = new OpenFileDialog();
            if (openFileDialog.ShowDialog() == System.Windows.Forms.DialogResult.OK)
            {
                try
                {
                    Mat SourceMat = new Mat(openFileDialog.FileName, ImreadModes.AnyColor);
                    Image<Rgb, byte> source = new Image<Rgb, byte>(SourceMat.Bitmap);
                    oriantation = new OriantatioOnMap(source,parametrs, Compression, Diameter);
                    SetImage(source);
                }
                catch { }
            }
        }

        private void Size_TextChanged(object sender, System.Windows.Controls.TextChangedEventArgs e)
        {
            int ColumnCount;
            int RowCount;
            try
            {
                ColumnCount = Math.Min(220, Convert.ToInt32(ColumnCountTextBox.Text));
                RowCount = Math.Min(220, Convert.ToInt32(RowCountTextBox.Text));
            }
            catch
            {
                return;
            }
            GridSize = new System.Drawing.Point(ColumnCount, RowCount);
            if (oriantation != null)
                SetImage(oriantation.Map);
        }

        private void CompressionTextBox_TextChanged(object sender, System.Windows.Controls.TextChangedEventArgs e)
        {            
            try
            {
                Compression = Math.Min(220, Convert.ToInt32(CompressionTextBox.Text));
            }
            catch
            {
                return;
            }
        }

        private void DiameterTextBox_TextChanged(object sender, System.Windows.Controls.TextChangedEventArgs e)
        {
            try
            {
                Diameter = Math.Min(220, Convert.ToInt32(DiameterTextBox.Text));
            }
            catch
            {
                return;
            }
        }

        private void StartFlightButton_Click(object sender, RoutedEventArgs e)
        {
            aircraft = new Aircraft(WayAircraft, 0.05);
        }

        private void RemoveWayButton_Click(object sender, RoutedEventArgs e)
        {
            WayAircraft.Clear();
            if (imagefromAircraft != null)
                SetImageAircraft(imagefromAircraft);
        }

        private void RemovelastPointButton_Click(object sender, RoutedEventArgs e)
        {
            WayAircraft.RemoveAt(WayAircraft.Count - 1);
            if (imagefromAircraft != null)
                SetImageAircraft(imagefromAircraft);
        }
        private void NextStepButton_Click(object sender, RoutedEventArgs e)
        {
            int width, height,k;
            double uniquenessThreshold;
            try
            {
                width = Convert.ToInt32(CameraWidthTextBox.Text);
                height = Convert.ToInt32(CameraHeightTextBox.Text);
                k = Convert.ToInt32(kTextBox.Text);
                uniquenessThreshold = Convert.ToDouble(UniquenessThresholdTextBox.Text);
            }
            catch
            { return; }
            Image<Rgb, byte> SubMap = 
                aircraft.NextLocate(imagefromAircraft/*oriantation.Map*/, new System.Drawing.Point(width, height));
            SetImage(oriantation.ShowMatches(SubMap, k, uniquenessThreshold,parametrs));

        }
    }
}
