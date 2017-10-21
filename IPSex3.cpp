#include <stdio.h>
#include <ctime>
#include <cilk/cilk.h>
#include <cilk/reducer_opadd.h>
#include <chrono>

using namespace std::chrono;

// количество строк в исходной квадратной матрице
const int MATRIX_SIZE = 1500;

/// Функция InitMatrix() заполняет переданную в качестве 
/// параметра квадратную матрицу случайными значениями
/// matrix - исходная матрица СЛАУ
void InitMatrix(double** matrix)
{
	for (int i = 0; i < MATRIX_SIZE; ++i)
	{
		matrix[i] = new double[MATRIX_SIZE + 1];
	}

	for (int i = 0; i < MATRIX_SIZE; ++i)
	{
		for (int j = 0; j <= MATRIX_SIZE; ++j)
		{
			matrix[i][j] = rand() % 10 + 1;//rand() % 2500 + 1;
		}
	}
}

/// Функция SerialGaussMethod() решает СЛАУ методом Гаусса 
/// matrix - исходная матрица коэффиициентов уравнений, входящих в СЛАУ,
/// последний столбей матрицы - значения правых частей уравнений
/// rows - количество строк в исходной матрице
/// result - массив ответов СЛАУ
double SerialGaussMethod(double **matrix, const int rows, double* result, double durS)
{
	int k;
	double koef;
	duration<double> duration; /// Перемнная для измерения времени

							   //Замеряем время для прямого хода метода Гаусса
	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	// прямой ход метода Гаусса
	for (k = 0; k < rows; ++k)
	{
		//
		for (int i = k + 1; i < rows; ++i)
		{
			koef = -matrix[i][k] / matrix[k][k];

			for (int j = k; j <= rows; ++j)
			{
				matrix[i][j] += koef * matrix[k][j];
			}
		}
	}
	high_resolution_clock::time_point t2 = high_resolution_clock::now();
	duration = (t2 - t1);
	durS = duration.count();
	printf("Duration is: %lf seconds (Serial code)\n", durS); // Выводим время работы прямого хода

															  // обратный ход метода Гаусса
	result[rows - 1] = matrix[rows - 1][rows] / matrix[rows - 1][rows - 1];

	for (k = rows - 2; k >= 0; --k)
	{
		result[k] = matrix[k][rows];

		//
		for (int j = k + 1; j < rows; ++j)
		{
			result[k] -= matrix[k][j] * result[j];
		}

		result[k] /= matrix[k][k];
	}
	return durS;
}
		/// Функция ParallelGaussMethod() решает СЛАУ методом Гаусса 
double ParallelGaussMethod(double **matrix, const int rows, double* result, double durP)
{
	int k;
	double koef;
	duration<double> duration; // Перемнная для измерения времени
					
	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	// прямой ход метода Гаусса
	for (k = 0; k < rows; ++k)
	{
		//использование cilk_for во внутреннем цикле прямого хода
		for (int i = k + 1; i < rows; ++i)
		{
			koef = -matrix[i][k] / matrix[k][k];
			cilk_for(int j = k; j <= rows; ++j)
			{
				matrix[i][j] += koef * matrix[k][j];
			}
		}
	}
	high_resolution_clock::time_point t2 = high_resolution_clock::now();
	duration = (t2 - t1);
	printf("Duration is: %lf seconds\n", duration.count()); // Выводим время работы прямого хода

	// обратный ход метода Гаусса
	result[rows - 1] = matrix[rows - 1][rows] / matrix[rows - 1][rows - 1];

	for (k = rows - 2; k >= 0; --k)
	{
		result[k] = matrix[k][rows];
		//
		cilk_for(int j = k + 1; j < rows; ++j)
			result[k] -= matrix[k][j] * result[j];


		result[k] /= matrix[k][k];
	}
}//
		for (int i = k + 1; i < rows; ++i)
		{
			koef = -matrix[i][k] / matrix[k][k];

			
			for (int j = k; j <= rows; ++j)
			{
				matrix[i][j] += koef * matrix[k][j];
			}
		}
	}
	high_resolution_clock::time_point t2 = high_resolution_clock::now();
	duration = (t2 - t1);
	printf("Duration is: %lf seconds\n", duration.count());
	// обратный ход метода Гаусса
	result[rows - 1] = matrix[rows - 1][rows] / matrix[rows - 1][rows - 1];

	for (k = rows - 2; k >= 0; --k)
	{
		result[k] = matrix[k][rows];
		cilk::reducer_opadd <double> res(result[k]);
		//
		for (int j = k + 1; j < rows; ++j)
		{
			res -= matrix[k][j] * result[j];
		}
		result[k] = res.get_value();
		result[k] /= matrix[k][k];
	}
	return durP;
}


int main()
{
	srand((unsigned)time(0));


	int i;
	double dS, dP;
	// кол-во строк в матрице, приводимой в качестве примера
	//const int test_matrix_lines = 4;//Изменяем размерность на MATRIX_SIZE

	double **test_matrix = new double*[MATRIX_SIZE];

	// цикл по строкам
	for (i = 0; i < MATRIX_SIZE; ++i)
	{
		// (test_matrix_lines + 1)- количество столбцов в тестовой матрице,
		// последний столбец матрицы отведен под правые части уравнений, входящих в СЛАУ
		test_matrix[i] = new double[MATRIX_SIZE + 1];
	}

	// массив решений СЛАУ
	double *result = new double[MATRIX_SIZE];
	double *rP = new double[MATRIX_SIZE];
	// инициализация тестовой матрицы
	/*test_matrix[0][0] = 2; test_matrix[0][1] = 5;  test_matrix[0][2] = 4;  test_matrix[0][3] = 1;  test_matrix[0][4] = 20;
	test_matrix[1][0] = 1; test_matrix[1][1] = 3;  test_matrix[1][2] = 2;  test_matrix[1][3] = 1;  test_matrix[1][4] = 11;
	test_matrix[2][0] = 2; test_matrix[2][1] = 10; test_matrix[2][2] = 9;  test_matrix[2][3] = 7;  test_matrix[2][4] = 40;
	test_matrix[3][0] = 3; test_matrix[3][1] = 8;  test_matrix[3][2] = 9;  test_matrix[3][3] = 2;  test_matrix[3][4] = 37;*/
	InitMatrix(test_matrix);

	dS = SerialGaussMethod(test_matrix, MATRIX_SIZE, result, dS);
	dP = ParallelGaussMethod(test_matrix, MATRIX_SIZE, rP, dP);
	for (i = 0; i < MATRIX_SIZE; ++i)
	{
		delete[]test_matrix[i];
	}
	printf("Acceleration = %lf\n", dS / dP);
	printf("Solution\n");
	for (i = 0; i < MATRIX_SIZE; ++i)
	{
		printf("x(%d) = %lf              x(%d) = %lf \n", i, result[i], i, resultP[i]);
	}

	delete[] result;
	delete[] rP;
	return 0;
}
