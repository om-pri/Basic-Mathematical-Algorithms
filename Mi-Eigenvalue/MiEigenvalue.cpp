#include <iostream>
#include <iomanip>
#include "math.h"

using namespace std;

//ȫ�ֱ���
int ii;	//��������������
double a[501], b = 0.16, c = -0.064;	//����A
double Ua[501], Ub[500], Uc[499], Lb[500], Lc[499];	//LU�ֽ���

/*�ݷ���KΪ������������pΪƽ����*/
double mi(int K, double p)
{
	//��������
	double u[501], u1[501], y[501];
	double beta, last_beta = 0, D, d;

	//����ƽ����p�������ĶԽ���Ԫ��
	for (int i = 1; i <= 501; i++)
		a[i - 1] = (1.64 - 0.024*i)*sin(0.2*i) - 0.64*exp(0.1 / double(i)) - p;

	//��ʼ����
	for (int i = 0; i < 501; i++)
		u[i] = 1;

	for (ii = 0; ii < K; ii++)
	{
		D = 0;
		d = 0;
		beta = 0;

		for (int i = 0; i < 501; i++)
			D += u[i] * u[i];
		
		//��һ������u
		d = sqrt(D);
		for (int i = 0; i < 501; i++)
			y[i] = u[i] / d;

		for (int i = 0; i < 501; i++)
		{
			if (i == 0)
				u1[i] = a[i] * y[i] + b*y[i + 1] + c*y[i + 2];
			else if (i == 1)
				u1[i] = b*y[i - 1] + a[i] * y[i] + b*y[i + 1] + c*y[i + 2];
			else if (i == 499)
				u1[i] = c*y[i - 2] + b*y[i - 1] + a[i] * y[i] + b*y[i + 1];
			else if (i == 500)
				u1[i] = c*y[i - 2] + b*y[i - 1] + a[i] * y[i];
			else
				u1[i] = c*y[i - 2] + b*y[i - 1] + a[i] * y[i] + b*y[i + 1] + c*y[i + 2];
		}

		for (int i = 0; i < 501; i++)
			beta += y[i] * u1[i];

		if (abs(beta - last_beta) / abs(beta) <= 1e-12)	//�ж��Ƿ����������ֹ����
			break;

		for (int i = 0; i < 501; i++)
			u[i] = u1[i];

		last_beta = beta;
	}
	return beta;
}

/*LU�ֽⷨ,pΪƽ����*/
void LU(double p)
{
	for (int i = 1; i <= 501; i++)
		a[i - 1] = (1.64 - 0.024*i)*sin(0.2*i) - 0.64*exp(0.1 / double(i)) - p;

	for (int i = 0; i < 501; i++)
	{
		if (i == 0){
			Ua[i] = a[i];
			Ub[i] = b;
			Uc[i] = c;
		}
		else if (i == 1){
			Lb[i - 1] = b / Ua[i - 1];
			Ua[i] = a[i] - Lb[i - 1] * Ub[i - 1];
			Ub[i] = b - Lb[i - 1] * Uc[i - 1];
			Uc[i] = c;
		}
		else if (i == 499){
			Lc[i - 2] = c / Ua[i - 2];
			Lb[i - 1] = (b - Lc[i - 2] * Ub[i - 2]) / Ua[i - 1];
			Ua[i] = a[i] - Lc[i - 2] * Uc[i - 2] - Lb[i - 1] * Ub[i - 1];
			Ub[i] = b - Lb[i - 1] * Uc[i - 1];
		}
		else if (i == 500){
			Lc[i - 2] = c / Ua[i - 2];
			Lb[i - 1] = (b - Lc[i - 2] * Ub[i - 2]) / Ua[i - 1];
			Ua[i] = a[i] - Lc[i - 2] * Uc[i - 2] - Lb[i - 1] * Ub[i - 1];
		}
		else{
			Lc[i - 2] = c / Ua[i - 2];
			Lb[i - 1] = (b - Lc[i - 2] * Ub[i - 2]) / Ua[i - 1];
			Ua[i] = a[i] - Lc[i - 2] * Uc[i - 2] - Lb[i - 1] * Ub[i - 1];
			Ub[i] = b - Lb[i - 1] * Uc[i - 1];
			Uc[i] = c;
		}
	}
}

/*���ݷ���KΪ����������*/
double in_mi(int K)
{
	//��������
	double u[501], u1[501], y[501], x[501];
	double beta, last_beta = 0, D, d;

	//��ʼ����
	for (int i = 0; i < 501; i++)
		u[i] = 1;

	for (ii = 0; ii < K; ii++)
	{
		D = 0;
		d = 0;
		beta = 0;

		for (int i = 0; i < 501; i++)
			D += u[i] * u[i];

		//��һ������u
		d = sqrt(D);
		for (int i = 0; i < 501; i++)
			y[i] = u[i] / d;

		//�ش�����
		for (int i = 0; i < 501; i++)
		{
			if (i == 0)
				x[i] = y[i];
			else if (i == 1)
				x[i] = y[i] - Lb[i - 1] * x[i - 1];
			else
				x[i] = y[i] - Lc[i - 2] * x[i - 2] - Lb[i - 1] * x[i - 1];
		}
		for (int i = 500; i >=0; i--)
		{
			if (i == 500)
				u1[i] = x[i]/Ua[i];
			else if (i == 499)
				u1[i] = (x[i] - Ub[i] * u1[i + 1]) / Ua[i];
			else
				u1[i] = (x[i] - Ub[i] * u1[i + 1] - Uc[i] * u1[i + 2]) / Ua[i];
		}

		for (int i = 0; i < 501; i++)
			beta += y[i] * u1[i];

		if (abs(beta - last_beta) / abs(beta) <= 1e-12)	//�ж��Ƿ����������ֹ����
			break;

		for (int i = 0; i < 501; i++)
			u[i] = u1[i];

		last_beta = beta;
	}
	return beta;
}

void main()
{
	double lam, temp, lam1, lambda_1, lambda_501, temp1, lambda_s;
	double detA = 1, condA = 1, p = 0;
	double u_k[39], lambda_k[39], temp2;
	int record_ii, record1_ii;
	
	int K;
	cout << "�������������������������������Դ���ֹ�� K= ";
	cin >> K;

	//�ݷ�����ĩ��������ֵ
	lam = mi(K, p);
	record_ii = ii;
	
	p = lam;
	temp = mi(K, p);
	lam1 = temp + lam;
	record1_ii = ii;
	
	if (lam > lam1){
		lambda_1 = lam1;
		lambda_501 = lam;
		cout << " ��1Ϊ��" << setprecision(12) << lambda_1 << "  ��������Ϊ��" << record1_ii << endl;
		cout << " ��501Ϊ��" << setprecision(12) << lambda_501 << "  ��������Ϊ��" << record_ii << endl;
	}
	else{
		lambda_1 = lam;
		lambda_501 = lam1;
		cout << " ��1Ϊ��" << setprecision(12) << lambda_1 << "  ��������Ϊ��" << record_ii << endl;
		cout << " ��501Ϊ��" << setprecision(12) << lambda_501 << "  ��������Ϊ��" << record1_ii << endl;
	}

	//���ݷ���ģ��С����ֵ
	LU(p = 0);	//���Ƚ���LU�ֽ�
	temp1 = in_mi(K);
	lambda_s = 1.0 / temp1;
	cout << " ��sΪ��" << setprecision(12) << lambda_s << "  ��������Ϊ��" << ii << endl;

	//�ڶ��ʣ����ik
	for (int i = 0; i < 39; i++)
	{
		u_k[i] = lambda_1 + ((i + 1)*(lambda_501 - lambda_1)) / 40;
		LU(p = u_k[i]);
		temp2 = in_mi(K);
		lambda_k[i] = 1.0 / temp2 + u_k[i];
		cout << " ��i" << i + 1 << "Ϊ��" << setprecision(12) << lambda_k[i] << "  ��������Ϊ��" << ii << endl;
	}

	//�����ʣ���A���׷���������������ʽֵ
	condA = (abs(lambda_1)>abs(lambda_501) ? abs(lambda_1) : abs(lambda_501)) / abs(lambda_s);
	cout << " A�ģ��׷�����������Ϊ��" << setprecision(12) << condA << endl;
	LU(p = 0);
	for (int i = 0; i < 501; i++)
		detA = detA*Ua[i];
	cout << " ����A����ʽֵΪ��" << setprecision(12) << detA << endl;
}
