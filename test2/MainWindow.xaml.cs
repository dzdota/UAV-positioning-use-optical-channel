
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

using System.Collections.Generic;
using System.Windows.Forms;
using SD = System.Drawing;
using SW = System.Windows;

namespace UVAPositioning
{
    /// <summary>
    /// Interaction logic for MainWindow.xaml
    /// </summary>
    public partial class MainWindow : Window
    {
        List<System.Windows.Point> WayAircraft = new List<System.Windows.Point>();
        OriantatioOnMap oriantation;
        Aircraft aircraft = null;
        Navigation navigation = null;
        Image<Rgb, byte> imagefromAircraft = null;
        SIFTParametrs parametrs;
        SD.Point GridSize;
        double Compression;
        double Diameter;
        string aircraftIconPath = "aircaftIcon.png";
        Icon baseIcon;
        Icon lostGPSIcon;

        object temp;


        public MainWindow()
        {

            InitializeComponent();
            parametrs = new SIFTParametrs()
            {
                nFeatures = 0,
                nOctaveLayers = 5,
                contrastThreshold = 0.04,
                edgeThreshold = 10,
                sigma = 1.6
            };

            /*
            Bitmap awad = new Bitmap("gps-disconnected.png");
            for (int i = 0; i < awad.Width; i++)
                for (int j = 0; j < awad.Height; j++)
                {
                    var col = awad.GetPixel(i, j);
                    if (col.A > 200 ||( col.R == 0 && col.G == 0 && col.B == 0))
                        awad.SetPixel(i, j, Color.FromArgb(0, 0, 255));
                    else
                        awad.SetPixel(i, j, Color.FromArgb(0, 0, 0, 0));


                }
            awad.Save("lost-signalIcon.png", SD.Imaging.ImageFormat.Png);
            */
        }

        public void SetImageAircraft(Image<Rgb, byte> image)
        {
            /*double Transform = Math.Max(oriantation.Map.Width / imageBox1.Width,
                oriantation.Map.Height / imageBox1.Height);
            image = image.Resize(Transform, Inter.Area);*/
            image = ImageTransform.SetLine(image, WayAircraft, new MCvScalar(0, 0, 255));
            if(aircraft != null)
            {
                image = ImageTransform.SetAircraft(image, aircraft);
            }
            if (baseIcon != null)
            {
                image = ImageTransform.SetIcon(image, baseIcon.IconImage, ImageTransform.SWtoSD(baseIcon.Coordinate));
            }
            if (lostGPSIcon != null)
            {
                image = ImageTransform.SetIcon(image, lostGPSIcon.IconImage, ImageTransform.SWtoSD(lostGPSIcon.Coordinate));
            }
            if (baseIcon != null) { }
            var imgbrush = new BitmapImage();
            imgbrush.BeginInit();
            var ms = new MemoryStream();
            image.Bitmap.Save(ms, SD.Imaging.ImageFormat.Bmp);
            imgbrush.StreamSource = ms;
            imgbrush.CreateOptions = BitmapCreateOptions.PreservePixelFormat;
            imgbrush.EndInit();
            imageBoxfromAircraft.Source = imgbrush;
        }
        public void SetMap(Image<Rgb, byte> image)
        {
            if (ShowGridCheckBox.IsChecked == true)
                image = ImageTransform.SetGrid(image, GridSize, oriantation.Map.Size);
            if (ShowKeyPointCheckBox.IsChecked == true && oriantation != null)
                image = ImageTransform.SetKeyPont(image, oriantation, new Bgr(255, 0, 0));
            var imgbrush = new BitmapImage();
            imgbrush.BeginInit();
            var ms = new MemoryStream();
            image.Bitmap.Save(ms, SD.Imaging.ImageFormat.Bmp);
            imgbrush.StreamSource = ms;
            imgbrush.CreateOptions = BitmapCreateOptions.PreservePixelFormat;
            imgbrush.EndInit();
            imageBox1.Source = imgbrush;
            imageBoxMiniMap.Source = imgbrush;
            image.Save("outputimage.bmp");
        }

        [Obsolete("Use OriantationOnMap")]
        public Bitmap DrawSift(Image<Rgb, byte> modelimage, Image<Rgb, byte> observedimage)
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
                SD.Rectangle rect = new SD.Rectangle(SD.Point.Empty, modelimage.Size);
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
                SD.Point[] points = Array.ConvertAll<PointF, SD.Point>(pts, SD.Point.Round);
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
                SD.Point p =
                    new SD.Point((int)observedKeyPoints[a[0].QueryIdx].Point.X,
                        (int)observedKeyPoints[a[0].QueryIdx].Point.Y);
                matrix[(int)p.Y * ypart / observedimage.Bitmap.Height, 
                    (int)p.X * xpart / observedimage.Bitmap.Width]++;
            }
            for (int x = 0; x < xpart; x++)
                for (int y = 0; y < ypart; y++)
                {
                    if (matrix[y, x] >= n)
                    {
                        SD.Rectangle rec = new SD.Rectangle(
                            (int)x * observedimage.Bitmap.Width / xpart,
                            (int)y * observedimage.Bitmap.Height / ypart,
                            observedimage.Bitmap.Width / xpart,
                            observedimage.Bitmap.Height / ypart);
                        CvInvoke.Rectangle(result, rec, new MCvScalar(0 ,0 , 255), 2);
                    }
                }
            return result;
        }

        private void ImageBox1_MouseLeftButtonDown(object sender, MouseButtonEventArgs e)
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
                    SetMap(source);
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
            GridSize = new SD.Point(ColumnCount, RowCount);
            if (oriantation != null)
                SetMap(oriantation.Map);
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
            aircraft = new Aircraft(WayAircraft, 0.05, new Image<Rgb, byte>(aircraftIconPath));
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
            double persent;
            try
            {
                width = Convert.ToInt32(CameraWidthTextBox.Text);
                height = Convert.ToInt32(CameraHeightTextBox.Text);
                k = Convert.ToInt32(kTextBox.Text);
                uniquenessThreshold = Convert.ToDouble(UniquenessThresholdTextBox.Text);
                persent = Convert.ToDouble(PersentFromFind.Text);
            }
            catch
            { return; }
            System.Windows.Point? locationInMap;
            if (navigation != null)
            {
                SW.Point location = navigation.NextLocation();
                WayAircraft.Add(location);
                Image<Rgb, byte> SubMap =
                    navigation.GetPhoto(imagefromAircraft, location, new SW.Point(width, height));
                SetMap(oriantation.ShowMatches(SubMap, k, uniquenessThreshold, GridSize.X, GridSize.Y, persent, parametrs, out locationInMap));
                if (locationInMap != null)
                    navigation.Location = location;
            }
            else  if (aircraft != null && aircraft.Locate < 1)
            {
                aircraft.FindLocate();
                Image<Rgb, byte> SubMap =
                    aircraft.GetPhoto(imagefromAircraft, new SD.Point(width, height));
                SetMap(oriantation.ShowMatches(SubMap, k, uniquenessThreshold, GridSize.X, GridSize.Y, persent, parametrs, out locationInMap));
            }

            SetImageAircraft(imagefromAircraft);

        }

        private void MenuItem_Click_SmoothGaussian(object sender, RoutedEventArgs e)
        {
            imagefromAircraft = imagefromAircraft.SmoothGaussian(3, 3, 34.3, 45.3);
            SetImageAircraft(imagefromAircraft);
        }

        private void MenuItem_Click_Remove_Noise(object sender, RoutedEventArgs e)
        {
            CvInvoke.FastNlMeansDenoisingColored(imagefromAircraft, imagefromAircraft);
            SetImageAircraft(imagefromAircraft);

        }

        private void MenuItem_Click_SmoothBilatral(object sender, RoutedEventArgs e)
        {
            imagefromAircraft = imagefromAircraft.SmoothBilatral(5, 30, 50);
            SetImageAircraft(imagefromAircraft);
        }

        private void MenuItem_Click_SmoothMedian(object sender, RoutedEventArgs e)
        {
            imagefromAircraft = imagefromAircraft.SmoothMedian(5);
            SetImageAircraft(imagefromAircraft);
        }

        private void MenuItem_Click_Add_Noise(object sender, RoutedEventArgs e)
        {
            Image<Rgb, byte> noise = new Image<Rgb, byte>(imagefromAircraft.Bitmap);
            noise.SetRandNormal(new MCvScalar(0,0,0), new MCvScalar(30, 30, 30));
            imagefromAircraft = imagefromAircraft + noise;
            SetImageAircraft(imagefromAircraft);

        }

        private void SetBase_Click(object sender, RoutedEventArgs e)
        {
            SD.Point screenPoint = Control.MousePosition;
            System.Windows.Point position = imageBoxfromAircraft.PointFromScreen(new System.Windows.Point(screenPoint.X, screenPoint.Y));
            double Transform = Math.Max(imagefromAircraft.Width / imageBoxfromAircraft.Width,
                imagefromAircraft.Height / imageBoxfromAircraft.Height);
            position = new System.Windows.Point(position.X * Transform,
                position.Y * Transform);
            baseIcon = new UVAPositioning.Icon(new Image<Rgb, byte>("base.png"), position);
            SetImageAircraft(imagefromAircraft);
            SetNavigation();
        }

        private void SetGPSLost_Click(object sender, RoutedEventArgs e)
        {
            SD.Point screenPoint = Control.MousePosition;
            System.Windows.Point position = imageBoxfromAircraft.PointFromScreen(new System.Windows.Point(screenPoint.X, screenPoint.Y));
            double Transform = Math.Max(imagefromAircraft.Width / imageBoxfromAircraft.Width,
                imagefromAircraft.Height / imageBoxfromAircraft.Height);
            position = new System.Windows.Point(position.X * Transform,
                position.Y * Transform);
            lostGPSIcon = new UVAPositioning.Icon(new Image<Rgb, byte>("lost-signalIcon.png"),position);
            SetImageAircraft(imagefromAircraft);
            SetNavigation();
        }

        private void SetNavigation()
        {
            if (oriantation != null && baseIcon != null && lostGPSIcon != null) {

                System.Windows.Point cameraSize = new System.Windows.Point();
                try
                {
                    cameraSize.X = Convert.ToInt32(CameraWidthTextBox.Text);
                    cameraSize.Y = Convert.ToInt32(CameraHeightTextBox.Text);
                }
                catch
                { return; }
                navigation = new Navigation(oriantation, baseIcon, lostGPSIcon, cameraSize);
            }
        }

        private void ClearIcon_Click(object sender, RoutedEventArgs e)
        {
            lostGPSIcon = null;
            baseIcon = null;
            navigation = null;
            WayAircraft.Clear();
            SetImageAircraft(imagefromAircraft);
        }
    }
}
