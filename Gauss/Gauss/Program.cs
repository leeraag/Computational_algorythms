using System;
using System.Diagnostics;

namespace Gauss
{
    class Program
    {
        private const double e = 1e-12;                    // Погрешность вычислений
        static void Main(string[] args)
        {

            Console.Write("Введите количество уравнений: ");
            int n = Convert.ToInt32(Console.ReadLine());
            Console.Write("Введите начальное значение диапазона коэфициентов и свободных членов системы: ");
            int a = Convert.ToInt32(Console.ReadLine());
            Console.Write("Введите конечное значение диапазона коэфициентов и свободных членов системы: ");
            int b = Convert.ToInt32(Console.ReadLine());
            var matrixA = new double[n, n];             // Коэффициенты
            var matrixAStart = new double[n, n];             // Коэффициенты
            var matrixB = new double[n];                // Свободные члены
            Random rnd = new Random();
            // Заполнение коэффициентов
            for (int i = 0; i < n; ++i)
            {
                for (int j = 0; j < n; ++j)
                {
                    int t = rnd.Next(a, b);
                    matrixA[i, j] = t;
                    matrixAStart[i, j] = t;
                }
            }

            // Заполнение свободных членов
            for (int i = 0; i < n; ++i)
                matrixB[i] = rnd.Next(a, b);

            PrintSlau(matrixA, matrixB);

            bool flag = true;

            var mtrxX = Gauss(matrixA, matrixB, n, true);
            //PrintVector("Результат: ", mtrxX);
            // Уточнение результата
            while (flag)
            {
                for (int i = 0; i < n; ++i)
                    for (int j = 0; j < n; ++j)
                        matrixA[i, j] = matrixAStart[i, j];
                // вычисление вектора невязки
                var mtrxB = new double[n]; //  вектор B0
                for (int i = 0; i < n; ++i)
                {
                    mtrxB[i] = 0;
                    for (int j = 0; j < n; ++j)
                        mtrxB[i] += matrixA[i, j] * mtrxX[j]; // a*x0=b0
                }
                var mtrxR = new double[n]; // Вектор невязки
                for (int i = 0; i < n; ++i)
                {
                    mtrxR[i] = matrixB[i] - mtrxB[i];
                }
                // Вывод невязки
                PrintVector("Результат: ", mtrxX);
                PrintVector("Невязка: ", mtrxR);


                var mtrxDelX = new double[n]; // Вектор поправки

                mtrxDelX = Gauss(matrixA, mtrxR, n, true); //a*x0=r

                for (int i = 0; i < n; ++i)
                {
                    mtrxX[i] = mtrxX[i] + mtrxDelX[i]; // новое приближенное решение системы x1
                }
                // Проверка на то, стоит ли продолжать
                for (int i = 0; i < n; ++i)
                {
                    if (mtrxDelX[i] > e)
                    {
                        flag = true;
                        break;
                    }
                    flag = false;
                }
            }

            // Конечный результат
            PrintVector("Конечный результат: ", mtrxX);
            Comparison(a, b);
        }

        private static void PrintSlau(double[,] matrixA, double[] matrixB)
        {
            int n = matrixB.Length;
            //Console.WriteLine("Матрица");
            for (int i = 0; i < n; ++i)
            {
                for (int j = 0; j < n; ++j)
                {
                    Console.Write(Math.Round(matrixA[i, j], 4) + "*x" + j);
                    if (j < n - 1)
                        Console.Write("\t+ ");
                }
                Console.WriteLine("\t= " + matrixB[i]);
            }
        }

        private static void PrintVector(string text, double[] matrix)
        {
            Console.WriteLine(text);
            for (int i = 0; i < matrix.Length; ++i)
                Console.WriteLine(i + ": " + matrix[i]);
        }

        private static double[] Gauss(double[,] matrixA, double[] matrixB, int n, bool show)
        {
            var matrixX = new double[n];
            int k = 0; // Номер шага
            // Прямой ход – приведение СЛАУ к треугольному виду
            while (k < n)
            {
                if (show) Console.WriteLine("\n\n\t\t\tШаг " + k + "\n\n");
                //if (show) PrintSlau(matrixA, matrixB);
                // Деление k-того уравнения на ведущий коэфициент
                double temp = matrixA[k, k];
                if (Math.Abs(temp) == 0) continue; // для нулевого коэффициента пропустить (если элемент гл диагонали равен 0)
                for (int j = 0; j < n; j++)
                    matrixA[k, j] = matrixA[k, j] / temp;
                matrixB[k] = matrixB[k] / temp;

                // Исключение переменной xk. Вычитание из k+1-ого уравнения k-того уравнения, умноженного на ведущий коэфициент k+1-ого 
                for (int i = k; i < n; i++)
                {
                    double temp2 = matrixA[i, k];
                    if (i == k)
                    {
                        //if (show) Console.WriteLine("i==k");
                        //if (show) PrintSlau(matrixA, matrixB);
                        continue; // уравнение не вычитать само из себя
                    }
                    for (int j = 0; j < n; j++)
                        matrixA[i, j] = matrixA[i, j] - matrixA[k, j] * temp2;
                    matrixB[i] = matrixB[i] - matrixB[k] * temp2;
                    //if (show) PrintSlau(matrixA, matrixB);
                }
                if (show) PrintSlau(matrixA, matrixB);
                k++;
            }

            // Обратный ход Гаусса
            for (k = n - 1; k >= 0; k--)
            {
                double s = 0;
                for (int j = n - 1; j > k; j--)
                {
                    s += matrixX[j] * matrixA[k, j];
                }
                matrixX[k] = matrixB[k] - s;
            }

            return matrixX;
        }

        private static double[] GaussСhoice(double[,] matrixA, double[] matrixB, int n)
        {
            var matrixX = new double[n];
            int k = 0; // Номер шага
            int index;
            double max;
            // Прямой ход – приведение СЛАУ к треугольному виду
            while (k < n)
            {
                // Поиск строки с максимальным a[i][k]
                max = Math.Abs(matrixA[k, k]);
                index = k;
                for (int i = k + 1; i < n; i++)
                {
                    if (Math.Abs(matrixA[i, k]) > max)
                    {
                        max = Math.Abs(matrixA[i, k]);
                        index = i;
                    }
                }
                // Перестановка строк
                if (max == 0)
                {
                    // нет ненулевых диагональных элементов
                    Console.WriteLine("Решение получить невозможно из-за нулевого столбца " + index + " матрицы A");
                    return new double[] { 0 };
                }
                for (int j = 0; j < n; j++)
                {
                    (matrixA[k, j], matrixA[index, j]) = (matrixA[index, j], matrixA[k, j]);
                }
                (matrixB[k], matrixB[index]) = (matrixB[index], matrixB[k]);

                // Деление k-того уравнения на ведущий коэфициент
                double temp = matrixA[k, k];
                if (Math.Abs(temp) == 0) continue; // для нулевого коэффициента пропустить
                for (int j = 0; j < n; j++)
                    matrixA[k, j] = matrixA[k, j] / temp;
                matrixB[k] = matrixB[k] / temp;

                // Исключение переменной xk. Вычитание из k+1-ого уравнения k-того уравнения, умноженного на ведущий коэфициент k+1-ого 
                for (int i = k; i < n; i++)
                {
                    double temp2 = matrixA[i, k];
                    if (i == k) continue; // уравнение не вычитать само из себя

                    for (int j = 0; j < n; j++)
                        matrixA[i, j] = matrixA[i, j] - matrixA[k, j] * temp2;
                    matrixB[i] = matrixB[i] - matrixB[k] * temp2;
                }
                k++;
            }

            // Обратный ход Гаусса
            for (k = n - 1; k >= 0; k--)
            {
                double s = 0;
                for (int j = n - 1; j > k; j--)
                {
                    s += matrixX[j] * matrixA[k, j];
                }
                matrixX[k] = (matrixB[k] - s);
            }

            return matrixX;
        }

        private static void Comparison(int a, int b)
        {
            Console.WriteLine(
                "\n\n\tСравнение временных характеристик решения СЛАУ методом Гаусса и Методом Гаусса с выбором главного элемента\n");
            Console.WriteLine("Порядок\t\tГаусс\t\tГаусс с выбором");
            for (int i = 100; i <= 200; i += 10)
            {
                Stopwatch GaussTime = new Stopwatch(), GaussChoiceTime = new Stopwatch();

                var matrixA = new double[i, i];             // Коэффициенты
                var matrixB = new double[i];                // Свободные члены
                Random rnd = new Random();
                // Заполнение коэффициентов
                for (int k = 0; k < i; ++k)
                {
                    for (int j = 0; j < i; ++j)
                    {
                        matrixA[k, j] = rnd.Next(a, b);
                    }
                }
                // Заполнение свободных членов
                for (int k = 0; k < i; ++k)
                    matrixB[k] = rnd.Next(a, b);

                GaussTime.Start();
                Gauss(matrixA, matrixB, i, false);
                GaussTime.Stop();

                GaussChoiceTime.Start();
                GaussСhoice(matrixA, matrixB, i);
                GaussChoiceTime.Stop();

                Console.WriteLine(i + "\t\t" + GaussTime.ElapsedTicks + "\t\t" + GaussChoiceTime.ElapsedTicks);
            }
        }
    }
}