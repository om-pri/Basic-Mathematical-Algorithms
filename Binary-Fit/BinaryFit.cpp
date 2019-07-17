#include <iostream>
#include <cmath>

using namespace std;

#define k_max 10	//���ʱk�����ֵ
const double vec_t[6] = { 0, 0.2, 0.4, 0.6, 0.8, 1.0 };
const double vec_u[6] = { 0, 0.4, 0.8, 1.2, 1.6, 2.0 };
const double table[6][6] = { { -0.5, -0.34, 0.14, 0.94, 2.06, 3.5 },
							 { -0.42, -0.5, -0.26, 0.3, 1.18, 2.38 },
							 { -0.18, -0.5, -0.5, -0.18, 0.46, 1.42 },
							 { 0.22, -0.34, -0.58, -0.5, -0.1, 0.62 },
							 { 0.78, -0.02, -0.5, -0.66, -0.5, -0.02 },
							 { 1.5, 0.46, -0.26, -0.66, -0.74, -0.5 }
							};
/* ��ֵ����� */
class INTERPOLATION
{
public:
	INTERPOLATION(double, double, int);	//���캯��
	~INTERPOLATION(){};					//��������

	//���ݳ�Ա
	int N;								//Newton������������
	int k_min;							//�ﵽsigma����ʱ��С��kֵ
	double det_x[4];					//Newton���еĦ�x
	double C[k_max + 1][k_max + 1];		//�������ϵ������C
	double U[11][21];					//f(x,y)���ݱ�
	double error;						//ţ�ٷ�����ˮƽ
	double sigma;						//��С���������Ͼ���ˮƽ
	double sum_error;					//��С������������


	//��Ա����	
	void NewtonMethod(double x, double y, double &u, double &t);	//Newton���������Է�����
	void GaussMethod(double b[4], double a[4][4], double(&x)[4]);	//Gauss��Ԫ�������Է�����
	double DualInterpolation(double x, double y);					//��Ƭ˫���β�ֵ
	double SurfaceFitting(int k, double x, double y);				//�������
	void GaussInverse(double(&matrix)[k_max + 1][k_max + 1], int k);//Gauss��ȥ���������	
	double max(double x, double y);									//��������ֵȡ��
};

/* ���캯�� */
INTERPOLATION::INTERPOLATION(double _error, double _sigma, int _N)
{
	error = _error;
	sigma = _sigma;
	N = _N;
}

/* Newton���������Է����� */
void INTERPOLATION::NewtonMethod(double x, double y, double &u, double &t)
{
	double v = 1.0, w = 1.0;
	u = 1.0;
	t = 1.0;
	double f[4], ff[4][4];
	int k = 0;
	for (k = 0; k <= N; k++)
	{
		f[0] = -1.0*(0.5*cos(t) + u + v + w - x - 2.67);
		f[1] = -1.0*(t + 0.5*sin(u) + v + w - y - 1.07);
		f[2] = -1.0*(0.5*t + u + cos(v) + w - x - 3.74);
		f[3] = -1.0*(t + 0.5*u + v + sin(w) - y - 0.79);
		ff[0][0] = -0.5*sin(t);	 ff[0][1] = 1.0;		 ff[0][2] = 1.0;	  ff[0][3] = 1.0;
		ff[1][0] = 1.0;			 ff[1][1] = 0.5*cos(u);	 ff[1][2] = 1.0;	  ff[1][3] = 1.0;
		ff[2][0] = 0.5;			 ff[2][1] = 1.0;		 ff[2][2] = -sin(v);  ff[2][3] = 1.0;
		ff[3][0] = 1.0;			 ff[3][1] = 0.5;		 ff[3][2] = 1.0;	  ff[3][3] = cos(w);
		//Gauss��ȥ��
		GaussMethod(f, ff, det_x);
		if (max(det_x[3], max(det_x[2], max(det_x[0], det_x[1]))) / max(w, max(v, max(u, t))) <= error)
			break;
		else
		{
			t += det_x[0];
			u += det_x[1];
			v += det_x[2];
			w += det_x[3];
			if (k == N)
				printf("���ʧ��!\n");
		}
	}
}

/* Gauss��ȥ�������Է����� */
void INTERPOLATION::GaussMethod(double b[4], double a[4][4], double(&x)[4])
{
	int i, j, i_k;
	double tmp, m_ik;
	for (int k = 0; k<3; k++)
	{
		i_k = k;
		//ѡ�к�i_k
		for (i = k; i < 4; i++)
		{
			if (fabs(a[i_k][k]) < fabs(a[i][k]))
				i_k = i;
		}
		//��������Ԫ��
		for (j = k; j<4; j++)
		{
			tmp = a[k][j];
			a[k][j] = a[i_k][j];
			a[i_k][j] = tmp;
		}
		tmp = b[k];
		b[k] = b[i_k];
		b[i_k] = tmp;

		for (i = k + 1; i<4; i++)
		{
			m_ik = a[i][k] / a[k][k];
			for (j = k; j<4; j++)
				a[i][j] = a[i][j] - m_ik*a[k][j];
			b[i] -= m_ik*b[k];
		}
	}
	//�ش�����
	x[3] = b[3] / a[3][3];
	for (int k = 2; k >= 0; k--)
	{
		tmp = 0;
		for (j = k + 1; j<4; j++)
			tmp += a[k][j] * x[j];
		x[k] = (b[k] - tmp) / a[k][k];
	}
}

/* ��������ֵȡ�� */
double INTERPOLATION::max(double x, double y)
{
	return fabs(x)>fabs(y) ? fabs(x) : fabs(y);
}

/* ��Ƭ˫���β�ֵ */
double INTERPOLATION::DualInterpolation(double x, double y)
{
	int i, j, k, r, t1, t2;
	double z = 0;
	i = int(fabs((x / 0.2) + 0.5));
	j = int(fabs((y / 0.4) + 0.5));
	if (i == 0)  i = 1;
	if (i == 5)  i = 4;
	if (j == 0)  j = 1;
	if (j == 5)  j = 4;

	for (k = i - 1; k <= i + 1; k++)
	{
		for (r = j - 1; r <= j + 1; r++)
		{
			double sum = 1.0;
			sum *= table[k][r];
			for (t1 = i - 1; t1 <= i + 1; t1++)
			{
				if (t1 != k)
					sum *= (x - vec_t[t1]) / (vec_t[k] - vec_t[t1]);
			}
			for (t2 = j - 1; t2 <= j + 1; t2++)
			{
				if (t2 != r)
					sum *= (y - vec_u[t2]) / (vec_u[r] - vec_u[t2]);
			}
			z += sum;
		}
	}
	return z;
}

/* ������� */
double INTERPOLATION::SurfaceFitting(int k, double x, double y)
{
	double B[11][k_max + 1] = { 0 }, BT[k_max + 1][11] = { 0 };
	double G[21][k_max + 1] = { 0 }, GT[k_max + 1][21] = { 0 };
	double BTB[k_max + 1][k_max + 1] = { 0 }, GTG[k_max + 1][k_max + 1] = { 0 };
	double matrix1[k_max + 1][11] = { 0 }, matrix2[k_max + 1][21] = { 0 }, matrix3[21][k_max + 1] = { 0 };
	double p, sum;
	int i, j, t, r, s;

	for (i = 0; i <= 10; i++)
	{
		for (j = 0; j <= k; j++)
			B[i][j] = pow(0.08*i, j);
	}
	for (i = 0; i <= 20; i++)
	{
		for (j = 0; j <= k; j++)
			G[i][j] = pow(0.5 + 0.05*i, j);
	}
	for (i = 0; i <= 10; i++)
	{
		for (j = 0; j <= k; j++)
			BT[j][i] = B[i][j];
	}
	for (i = 0; i <= 20; i++)
	{
		for (j = 0; j <= k; j++)
			GT[j][i] = G[i][j];
	}

	for (i = 0; i <= k; i++)
	{
		for (j = 0; j <= k; j++)
		{
			for (sum = 0, t = 0; t <= 10; t++)
				sum += BT[i][t] * B[t][j];
			BTB[i][j] = sum;
		}
	}
	GaussInverse(BTB, k);	//��BTB�������
	for (i = 0; i <= k; i++)
	{
		for (j = 0; j <= k; j++)
		{
			for (sum = 0, t = 0; t <= 20; t++)
				sum += GT[i][t] * G[t][j];
			GTG[i][j] = sum;
		}
	}
	GaussInverse(GTG, k);	//��GTG�������
	for (i = 0; i <= k; i++)
	{
		for (j = 0; j <= 10; j++)
		{
			for (sum = 0, t = 0; t <= k; t++)
				sum += BTB[i][t] * BT[t][j];
			matrix1[i][j] = sum;
		}
	}
	for (i = 0; i <= 20; i++)
	{
		for (j = 0; j <= k; j++)
		{
			for (sum = 0, t = 0; t <= k; t++)
				sum += G[i][t] * GTG[t][j];
			matrix3[i][j] = sum;
		}
	}
	for (i = 0; i <= k; i++)
	{
		for (j = 0; j <= 20; j++)
		{
			for (sum = 0, t = 0; t <= 10; t++)
				sum += matrix1[i][t] * U[t][j];
			matrix2[i][j] = sum;
		}
	}
	for (i = 0; i <= k; i++)
	{
		for (j = 0; j <= k; j++)
		{
			for (sum = 0, t = 0; t <= 20; t++)
				sum += matrix2[i][t] * matrix3[t][j];
			C[i][j] = sum;
		}
	}
	p = 0;
	for (r = 0; r <= k; r++)
	{
		for (s = 0; s <= k; s++)
			p += C[r][s] * pow(x, r)*pow(y, s);
	}
	return p;
}

/* Gauss��ȥ��������� */
void INTERPOLATION::GaussInverse(double(&matrix)[k_max + 1][k_max + 1], int k)
{
	double matrixT[k_max + 1][k_max + 1] = { 0 };
	double b[k_max + 1][k_max + 1] = { 0 }, tmp, m_ik;
	int i, j, t, p, i_k;
	for (i = 0; i <= k_max; i++)
	{
		for (j = 0; j <= k_max; j++)
		{
			if (i == j)
				b[i][j] = 1.0;
			else
				b[i][j] = 0;
		}
	}

	for (t = 0; t<k; t++)
	{
		i_k = t;
		for (i = t + 1; i <= k; i++)
		{
			if (fabs(matrix[i_k][t]) < fabs(matrix[i][t]))
				i_k = i;
		}
		for (j = t; j <= k; j++)
		{
			tmp = matrix[i_k][j];
			matrix[i_k][j] = matrix[t][j];
			matrix[t][j] = tmp;
		}

		for (p = 0; p <= k; p++)
		{
			tmp = b[t][p];
			b[t][p] = b[i_k][p];
			b[i_k][p] = tmp;
		}

		for (i = t + 1; i <= k; i++)
		{
			m_ik = matrix[i][t] / matrix[t][t];
			for (j = t + 1; j <= k; j++)
				matrix[i][j] -= m_ik*matrix[t][j];
			for (p = 0; p <= k; p++)
				b[i][p] -= m_ik*b[t][p];
		}
	}

	for (p = 0; p <= k; p++)
	{
		matrixT[k][p] = b[k][p] / matrix[k][k];
		for (i = k - 1; i >= 0; i--)
		{
			for (tmp = 0, j = i + 1; j <= k; j++)
				tmp += matrix[i][j] * matrixT[j][p];
			matrixT[i][p] = (b[i][p] - tmp) / matrix[i][i];
		}
	}
	for (i = 0; i <= k; i++)
	{
		for (j = 0; j <= k; j++)
			matrix[i][j] = matrixT[i][j];
	}
}

/* ������ */
void main()
{
	double x, y, t, u;
	INTERPOLATION ipfit(1e-12, 1e-7, 100);	//���ʵ����
	cout << "(����xi,yi,f(xi,yi)) (i=0,1...10;j=0,1...20):" << endl;
	for (int i = 0; i <= 10; i++)
	{
		for (int j = 0; j <= 20; j++)
		{
			x = 0.08*i;
			y = 0.5 + 0.05*j;
			ipfit.NewtonMethod(x, y, u, t);					//Newton������
			ipfit.U[i][j] = ipfit.DualInterpolation(t, u);	//��Ƭ˫���β�ֵ
			printf("x%d=%f, y%d=%f, f(x%d, y%d)=%.12le\n", i, x, j, y, i, j, ipfit.U[i][j]);
		}
	}
	cout << "ѡ����̵�k,��ֵ�ֱ�Ϊ��" << endl;
	for (int k = 0; k <= k_max; k++)
	{
		ipfit.sum_error = 0;
		for (int i = 0; i <= 10; i++)				//�������
		{
			for (int j = 0; j <= 20; j++)
				ipfit.sum_error += pow(ipfit.SurfaceFitting(k, 0.08*i, 0.5 + 0.05*j) - ipfit.U[i][j], 2);
		}
		printf("%d %.12le\n", k, ipfit.sum_error);
		if (ipfit.sum_error <= ipfit.sigma)			//�ж��Ƿ�ﵽ��Ͼ���Ҫ��
		{
			ipfit.k_min = k;
			printf("�ﵽ����Ҫ���k,��ֵ�ֱ�Ϊ��\nk=%d,��=%.12le\n", ipfit.k_min, ipfit.sum_error);
			cout << "ϵ��Crs(r=0,1,...,k;s=0,1,...,k)Ϊ:" << endl;
			for (int i = 0; i <= k; i++)
			{
				for (int j = 0; j <= k; j++)
					printf("Crs[%d][%d]=%.12le\n", i, j, ipfit.C[i][j]);
			}
			break;
		}
		else if (k == k_max)
			cout << "���ʧ�ܣ�δ�ﵽ����Ҫ��!" << endl;
	}
	//�۲�p(x,y)�ƽ�f(x,y)��Ч��
	cout << "����(x*[i],y*[j],f(x*[i],y*[j]),p(x*[i],y*[j]))(i=1,2,...,8;j=1,2,...,5):" << endl;
	for (int i = 1; i <= 8; i++)
	{
		for (int j = 1; j <= 5; j++)
		{
			x = 0.1*i;
			y = 0.5 + 0.2*j;
			ipfit.NewtonMethod(x, y, u, t);
			printf("x*[%d]=%f,y*[%d]=%f\n", i, x, j, y);
			printf("f(x*[%d],y*[%d])=%.12le, ", i, j, ipfit.DualInterpolation(t, u));
			printf("p(x*[%d],y*[%d])=%.12le\n", i, j, ipfit.SurfaceFitting(ipfit.k_min, x, y));
		}
	}
}
