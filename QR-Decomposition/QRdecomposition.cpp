#include <iostream>
#include <iomanip>
#include "math.h"

using namespace std;

//宏定义及全局变量
#define N 10
#define ERROR 1e-12
double a[N][N];
int m = N;
double lambda[N][2];	//矩阵A的全部特征值

/*初始化矩阵A*/
void init()
{
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			if (i == j)
				a[i][j] = 1.52*cos((i + 1) + 1.2*(j + 1));
			else
				a[i][j] = sin(0.5*(i + 1) + 0.2*(j + 1));
		}
	}
}

/*矩阵的拟上三角化*/
void Hessenberg()
{
	int count, r;
	double d, c, h, t;
	double u[N], p[N], q[N], w[N];

	for (r = 0; r < N-2; r++)
	{
		count = 0;
		for (int i = r + 2; i < N; i++)
		{
			if (fabs(a[i][r])>ERROR)
				count++;
		}
		if (count == 0)
			continue;
		else
		{
			d = 0;
			for (int i = r + 1; i < N; i++)
				d += a[i][r] * a[i][r];
			d = sqrt(d);
			c = (a[r + 1][r] > 0) ? (-d) : d;
			h = c*c - c*a[r + 1][r];

			for (int i = 0; i < N; i++)
			{
				if (i < r + 1)
					u[i] = 0;
				else if (i == r + 1)
					u[i] = a[i][r] - c;
				else
					u[i] = a[i][r];
			}

			for (int i = 0; i < N; i++)
			{
				p[i] = 0;
				for (int j = 0; j < N; j++)
					p[i] += a[j][i] * u[j];
				p[i] = p[i] / h;
			}
			for (int i = 0; i < N; i++)
			{
				q[i] = 0;
				for (int j = 0; j < N; j++)
					q[i] += a[i][j] * u[j];
				q[i] = q[i] / h;
			}

			t = 0;
			for (int i = 0; i < N; i++)
				t += p[i] * u[i];
			t = t / h;

			for (int i = 0; i < N; i++)
				w[i] = q[i] - t*u[i];

			for (int i = 0; i < N; i++)
				for (int j = 0; j < N; j++)
					a[i][j] = a[i][j] - w[i] * u[j] - u[i] * p[j];
		}
	}
}

/*矩阵M的QR分解与矩阵A的计算*/
void QR()
{
	int count;
	double s, det, d, c, h, t;
	double u[10], v[10], p[10], q[10], w[10], b[10][10] = { 0 };

	s = a[m - 2][m - 2] + a[m - 1][m - 1];
	det = a[m - 2][m - 2] * a[m - 1][m - 1] - a[m - 1][m - 2] * a[m - 2][m - 1];
	//计算出矩阵M
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < m; j++)
		{
			for (int r = 0; r < m; r++)
				b[i][j] += a[i][r] * a[r][j];
			b[i][j] = b[i][j] - s*a[i][j] + det*(i == j);
		}
	}
	
	for (int r = 0; r < m - 1; r++)
	{
		count = 0;
		for (int i = r + 1; i < m; i++)
		{
			if (fabs(b[i][r])>ERROR)
				count++;
		}
		if (count == 0)
			continue;
		else
		{
			d = 0;
			for (int i = r; i < m; i++)
				d += b[i][r] * b[i][r];
			d = sqrt(d);
			c = (b[r][r] > 0) ? (-d) : d;
			h = c*c - c*b[r][r];

			for (int i = 0; i < m; i++)
			{
				if (i < r)
					u[i] = 0;
				else if (i == r)
					u[i] = b[i][r] - c;
				else
					u[i] = b[i][r];
			}

			for (int i = 0; i < m; i++)
			{
				v[i] = 0;
				for (int j = 0; j < m; j++)
					v[i] += b[j][i] * u[j];
				v[i] = v[i] / h;
			}

			for (int i = 0; i < m; i++)
				for (int j = 0; j < m; j++)
					b[i][j] = b[i][j] - u[i] * v[j];

			for (int i = 0; i < m; i++)
			{
				p[i] = 0;
				for (int j = 0; j < m; j++)
					p[i] += a[j][i] * u[j];
				p[i] = p[i] / h;
			}
			for (int i = 0; i < m; i++)
			{
				q[i] = 0;
				for (int j = 0; j < m; j++)
					q[i] += a[i][j] * u[j];
				q[i] = q[i] / h;
			}

			t = 0;
			for (int i = 0; i < m; i++)
				t += p[i] * u[i];
			t = t / h;

			for (int i = 0; i < m; i++)
				w[i] = q[i] - t*u[i];

			for (int i = 0; i < m; i++)
			{
				for (int j = 0; j < m; j++)
					a[i][j] = a[i][j] - w[i] * u[j] - u[i] * p[j];
			}
		}
	}
}

/*解一元二次方程*/
void root()
{
	double delt, delt_1, b_1;
	b_1 = -(a[m - 2][m - 2] + a[m - 1][m - 1]);
	delt = b_1*b_1 - 4 * (a[m - 2][m - 2] * a[m - 1][m - 1] - a[m - 1][m - 2] * a[m - 2][m - 1]);
	delt_1 = sqrt(fabs(delt));
	if (delt >= 0)  //实特征值
	{
		lambda[m - 1][0] = (-b_1 - delt_1) / 2;
		lambda[m - 2][0] = (-b_1 + delt_1) / 2;
	}
	else			//复特征值
	{
		lambda[m - 1][0] = (-b_1) / 2;
		lambda[m - 1][1] = (-delt_1) / 2;
		lambda[m - 2][0] = (-b_1) / 2;
		lambda[m - 2][1] = delt_1 / 2;
	}
}

/*带双步位移的QR分解法*/
void doubleQR()
{
	int k = 1;
loop3:
	if (fabs(a[m - 1][m - 2]) < ERROR)
	{
		lambda[m - 1][0] = a[m - 1][m - 1];
		m = m - 1;
		goto loop4;
	}
	else
		goto loop5;
loop4:	
	if (m == 2)
	{
		root();
		return;	//计算成功
	}
	else if (m == 1)
	{
		lambda[m - 1][0] = a[m - 1][m - 1];
		return;	//计算成功
	}
	else
		goto loop3;
loop5:  
	if (fabs(a[m - 2][m - 3]) < ERROR)
	{
		root();
		m = m - 2;
		goto loop4;
	}
	else
		goto loop6;
loop6:
	if (k == 100) //最大迭代次数
	{
		cout << "计算失败，程序停止！" << endl;
		return;
	}
	else
		goto loop7;
loop7:
	QR();	//矩阵M的QR分解与矩阵A的计算
	k = k + 1;
	goto loop3;
}

/*列主元素Gauss消去法*/
void Gauss(double lambda_i)
{
	int i_k;
	double b, m, x[N] = {0};
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			if (i == j)
				a[i][j] = 1.52*cos((i + 1) + 1.2*(j + 1)) - lambda_i;
			else
				a[i][j] = sin(0.5*(i + 1) + 0.2*(j + 1));
		}
	}
	//消元过程
	for (int k = 0; k < N - 1; k++)
	{
		i_k = k;
		//选行号i_k
		for (int i = k+1; i < N; i++)
		{
			if (fabs(a[i][k])>fabs(a[i_k][k]))
				i_k = i;
		}
		//交换两行元素
		for (int j = k; j < N; j++)
		{
			b = a[k][j];
			a[k][j] = a[i_k][j];
			a[i_k][j] = b;
		}
		for (int i = k + 1; i < N; i++)
		{
			m = a[i][k] / a[k][k];
			for (int j = k + 1; j < N; j++)
				a[i][j] -= m*a[k][j];
		}
	}
	//回代过程
	x[9] = 1;
	for (int k = N - 2; k >= 0; k--)
	{
		for (int j = N - 1; j>k; j--)
			x[k] -= a[k][j] * x[j];
		x[k] /= a[k][k];
	}
	for (int i = 0; i < N; i++)
		cout << x[i] << endl;
}

void main()
{
	init();	  //初始化
	Hessenberg();	//拟上三角化
	cout << "拟上三角化矩阵A(n-1)为:" << endl;
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			if (fabs(a[i][j]) < ERROR)
				a[i][j] = 0;
			cout << setiosflags(ios::scientific) << setprecision(12) << a[i][j] << " ";
		}
		cout << endl;
	}
	doubleQR();	  //带双步位移的QR分解
	cout << "QR方法迭代结束后得到的矩阵为:" << endl;
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			if (fabs(a[i][j]) < ERROR)
				a[i][j] = 0;
			cout << setiosflags(ios::scientific) << setprecision(12) << a[i][j] << " ";
		}
		cout << endl;
	}
	cout << "矩阵的全部特征值为:" << endl;
	for (int i = 0; i < N; i++)
	{
		cout << setiosflags(ios::scientific) << setprecision(12) << "(" << lambda[i][0] << " ," << lambda[i][1] << ")" << endl;
	}

	for (int i = 0; i < N; i++)
	{
		if (lambda[i][1] == 0)
		{
			cout << setiosflags(ios::scientific) << setprecision(12) << "从属于" << lambda[i][0] << "的特征向量为：" << endl;
			Gauss(lambda[i][0]);
		}			
	}
}
