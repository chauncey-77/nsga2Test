#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <memory.h>
#include <vector>
#include <fstream>
#include <algorithm>
#include <cassert>
//#include "random.h"
#include "recombination.h"
using namespace std; //Ӧ����ȫ�֣����ù���ֱ����
double lowBound = 0;//�Ա���ȡֵ��Χ����
double uppBound = 1;//�Ա���ȡֵ��Χ����
int nvar = 30;//�Ա���x�ĸ���
int nobj = 2;//Ŀ�꺯��y�ĸ���
int pops = 100;//��Ⱥ����
int max_gen = 500;//���е������������Բ����500����Ϊ�������и�ֵ
int seed = 109;//�������������ֵ����������
long rnd_uni_init;//64λ���з������Σ��������Ϊ�˴���random.h�е�rnd_uni��������0��1֮��������
char filename[1024];
FILE *fp1, *fp2;
//ʵ������SBX����Ĳ���
int etax = 20;
int etam = 20;
double realx, realm; // SBX������� and ����ʽ����
double realb = 0.9; // ����Ⱥ��ѡ�񹩽����0
int gID;
void execute(char *algorithm);
void minfastsort(double x[], int idx[], int n, int m);
void objective(vector <double> &x_var, vector <double> &y_obj);
double rnd_uni(long *l);
struct point
{
	float x1, x2;
	float y1, y2;
};
bool cmp(point p1, point p2)
{
	if (p1.x1 < p2.x1) return true;
	return false;
}
class zc_NSGA2Ind
{
public:
	vector <double> x_var;//��̬���飬���ڴ���Ա���x
	vector <double> y_obj;//��̬���飬���ڴ��Ŀ��y
	int np;//��֧�������
	int sp;//֧����˵�����
	void rnd_init();//�����ʼ��
	void obj_eval();//����Ŀ�꺯��
	void show_objective();//�˺�������ʾĿ�꺯���������ڵ��ù�����û�õ������Ե�ʱ������������
	zc_NSGA2Ind();//���캯�������ô���ʱ�Զ����ô˺��������ڳ�ʼ��
	~zc_NSGA2Ind();//���������������ô˺����������������
	bool operator < (const zc_NSGA2Ind &ind2);//��"<"����,���¸��������壬ʹ"<"���ڸ���֮����бȽ�
	bool operator == (const zc_NSGA2Ind &ind2);//ͬ�ϣ���"=="���أ�
	void operator = (const zc_NSGA2Ind &ind2);//ͬ��
};

//���캯�������ô���ʱ�Զ����ô˺��������ڳ�ʼ��
zc_NSGA2Ind::zc_NSGA2Ind()
{
	for (int i = 0; i < nvar; i++)
	{
		x_var.push_back(0.0);
	}
	for (int j = 0; j < nobj; j++)
	{
		y_obj.push_back(0.0);
	}
}

//���������������ô˺����������������
zc_NSGA2Ind::~zc_NSGA2Ind()
{
	x_var.clear();
	y_obj.clear();
}

class zc_NSGA2
{
public:
	zc_NSGA2();//���캯�����ڱ�������ûʲô������
	virtual ~zc_NSGA2();//������������������ûʲô������
	void execute(int run);//ִ�к�����������ͨ���˺������ñ����еĸ�������
	void init_population();//��ʼ�����壬�������Ա�����Ŀ�꺯���ĳ�ʼ��
	int tour_selection();//���̶�ѡ�����
	void fill_union(zc_NSGA2Ind &ind);//��N������ŵ�offspring�У�����֤û���ظ�
	void rank_population();//��֧������
	void evaluation_density();//ӵ��������ѡ�񣬴˹����ڷ�֧������֮��
	void evolution();//���������Ӵ��Ĺ��̣���������ͱ���
	void save_front(char savefilename[1024]);//������Ŀ�꺯��ֵ�����POF�ļ�����
	void save_ps(char savefilename[1024]);//������Ŀ�꺯������Ӧ���Ա��������POS�ļ�����
	vector <zc_NSGA2Ind> population;//��Ⱥ����P���ϣ�ע��������ŵĸ���������CNSGA2Ind
	vector <zc_NSGA2Ind> offspring;//��ȺP���Ͻ�������Ӵ�Q���ϣ�ע��ͬ��
};

//���캯��
zc_NSGA2::zc_NSGA2() {}

//��������
zc_NSGA2::~zc_NSGA2() {}

//��"<"����,���¸��������壬ʹ"<"���ڸ���֮����бȽ�
bool zc_NSGA2Ind::operator < (const zc_NSGA2Ind &ind2)
{
	int flag2 = 0;
	for (int n = 0; n < nobj; n++)
	{
		if (ind2.y_obj[n] < y_obj[n])
			return false;
		if (ind2.y_obj[n] == y_obj[n])
			flag2++;
	}
	if (flag2 == nobj)
		return false;
	return true;
}

//��"=="���أ����¸�������˼��ʹ"=="���ڸ���֮����бȽ�
bool zc_NSGA2Ind::operator ==(const zc_NSGA2Ind &ind2)
{
	int flag = 0;
	for (int n = 0; n < nobj; n++)//nobj��Ŀ�꺯���ĸ���
	{
		if (ind2.y_obj[n] != y_obj[n])
			return false;
	}
	return true;
}

//��"="���أ����¸�������˼��ʹ"="�����ڸ���֮�临��
void zc_NSGA2Ind::operator = (const zc_NSGA2Ind &ind2)
{
	for (int n = 0; n <nvar; n++)
		x_var[n] = ind2.x_var[n];
	for (int n = 0; n < nobj; n++)
		y_obj[n] = ind2.y_obj[n];
	np = ind2.np;
}

int main()
{
	srand(time(NULL));
	execute("NSGA-II");
	return 0;
}

void execute(char *algorithm)
{
	int run = 1;
	zc_NSGA2 NSGA2;
	NSGA2.execute(run);
}

void zc_NSGA2::execute(int run)
{
	seed = (seed + 23) % 1377;//����һ���������
	rnd_uni_init = -(long)seed;//�����Ӳ���һ�������
	FILE *fp1, *fp2; /*�����ļ�ָ��*/
	FILE *fp;
	int gen = 1;//��һ��
	cout << "***********��" << gen << "��!**********" << endl;
	init_population();//���ó�ʼ������
	cout << endl;
	for (gen = 2; gen <= max_gen; gen++)
	{
		cout << "***********��" << gen << "��!**********" << endl;
		evolution();//�������������棬����
		rank_population();//����
		evaluation_density();//ӵ���ȱȽ�
		if (gen % max_gen == 0)
		{
			cout << "�쳿1";
			fopen_s(&fp1, "test1.txt", "w+"); /*����һ�������ļ�ֻд*/
			fopen_s(&fp2, "test2.txt", "w+"); /*����һ�������ļ�ֻд*/
			fopen_s(&fp, "data.txt", "w+");
			cout << "�쳿2";
			for (int n = 0; n < pops; n++)
			{
				for (int k = 0; k < nobj; k++)
				{
					fprintf(fp1, "%f ", population[n].y_obj[k]); /*�������ļ�дһ�ַ���*/
					cout << "POF = " << population[n].y_obj[k] << " ";
					cout << "�쳿3 : ";
				}
				fprintf(fp1, "\n");
			}
			for (int m = 0; m < pops; m++)
			{
				for (int s = 0; s < nvar; s++)
				{
					fprintf(fp2, "%f ", population[m].x_var[s]); /*�������ļ�дһ�ַ���*/
					cout << "POS = " << population[m].x_var[s] << " ";
					cout << "�쳿4 : ";
				}
				fprintf(fp2, "\n");
			}
			for (int n = 0; n < pops; ++n)
			{
				//for (int s = 0; s < nvar; ++s)
				//{
					//fprintf(fp, "%f ", population[n].x_var[s]);
				//}
				for (int s = 0; s < nobj; ++s)
				{
					fprintf(fp, "%f ", population[n].y_obj[s]);
				}
				fprintf(fp, "\n");
			}
			fclose(fp1); fclose(fp2);
			fclose(fp);

			/*ifstream in("data.txt");
			vector<point> ps;
			point p;
			while (in >> p.x1 >> p.x2 >> p.y1 >> p.y2)
			{
				ps.push_back(p);
			}
			sort(ps.begin(), ps.end(), cmp);
			in.close();
			ofstream out("data.txt");
			for (int i = 0; i < ps.size(); ++i)
			{
				out << ps[i].x1 << " " << ps[i].x2 << " " << ps[i].y1 << " " << ps[i].y2 << "\n";
			}
			out.close();*/
		}
	}
	population.clear();
	offspring.clear();
}

//��ʼ�����壬�������Ա�����Ŀ�꺯���ĳ�ʼ��
void zc_NSGA2::init_population()
{
	for (int n = 0; n < pops; n++)
	{
		zc_NSGA2Ind ind;//������һ�������
		ind.rnd_init();
		population.push_back(ind);
	}
}

//�����ʼ��
void zc_NSGA2Ind::rnd_init()
{
	for (int n = 0; n < nvar; n++)//nvarΪ�Ա����������˴�nvar=2
	{
		x_var[n] = lowBound + rnd_uni(&rnd_uni_init) * (uppBound - lowBound);
		printf("��%d�������ǣ�", n);
		std::cout << x_var[n] << " ";
		std::cout << "\n";
	}
	obj_eval();
	show_objective();
}

//����Ŀ�꺯��
void zc_NSGA2Ind::obj_eval()
{
	objective(x_var, y_obj);
}

//Ŀ�꺯��
void objective(vector <double> &x_var, vector <double> &y_obj)
{
	//y_obj[0] = (x_var[0] * x_var[0]) - x_var[1];
	//y_obj[1] = (x_var[0] + x_var[1] - 2)*x_var[0];

	y_obj[0] = x_var[0];
	double g = 0;
	for (int i = 1; i < nvar; ++i)
	{
		g = g + x_var[i];
	}
	g = g / (nvar - 1);
	g = 1 + 9 * g;
	double h = 1 - sqrt(x_var[0] / g);
	y_obj[1] = g * h;
}

//�˺�������ʾĿ�꺯��
void zc_NSGA2Ind::show_objective()
{
	for (int n = 0; n < nobj; n++)//nobj��Ŀ�꺯���ĸ���
	{
		std::cout << "Ŀ��: " << n << " ";
		std::cout << y_obj[n] << "\n";
	}
}

//���������Ӵ��Ĺ��̣���������ͱ���
void zc_NSGA2::evolution()
{
	for (int n = 0; n < pops; n++)
	{
		int p1, p2;
		p1 = tour_selection();//�����̶�ѡ��ķ���ֵ����P1
		while (1)//��ѭ����Ŀ������ѡ��������ͬ�ĸ��壬�Ա������н���
		{
			p2 = tour_selection();//�����̶�ѡ��ķ���ֵ����P2
			if (p1 != p2)
				break;
		}
		zc_NSGA2Ind child;
		//���ý��溯��(�⺯��)
		real_sbx_xover2(population[p1], population[p2], child);
		//���ñ��캯��(�⺯��)
		realmutation(child, 1.0 / nvar);
		//�����Ӵ�Ŀ��ֵ
		child.obj_eval();
		fill_union(child);
		fill_union(population[n]);
	}
}

int zc_NSGA2::tour_selection()
{
	int p1 = int(rnd_uni(&rnd_uni_init) * pops) % pops;//����һ�������
	int p2 = int(rnd_uni(&rnd_uni_init) * pops) % pops;//�����ڶ��������
	if (population[p1].np < population[p2].np)
		return p1;
	else
		return p2;
}

//��֧�����򣬼���offspring�����е�ÿһ������һ��rankֵ
void zc_NSGA2::rank_population()
{
	int size = offspring.size();//size = 2N
	int** cset;//����һ����ά����
	cset = new int*[size];
	for (int i = 0; i < size; i++)
		cset[i] = new int[size];
	//cset[2N][2N],cset[0][7]=96��ʾ��0����֧��ĵ�7������96
	int* rank = new int[size];
	//rank[2N]����ʾ��֧��ĵ㣬��rank[0]=2��ʾ��0���㱻������֧��
	cout << "size= " << size << "\n\n";
	for (int i = 0; i < size; i++)
	{
		rank[i] = 0;//֧���i����ĸ���Ϊ0
		offspring[i].np = -1;//֧���i����ĸ���Ϊ-1
		offspring[i].sp = 0;//��i������֧��ļ���Ϊ��
	}//��ʼ��
	for (int k = 0; k < size; k++)                           //����ò�ƿ��ԸĽ�һ��   xct
	{
		for (int j = 0; j < size; j++)
		{
			if (k != j)
			{
				if (offspring[j] < offspring[k] && !(offspring[j] == offspring[k]))//j֧��k
					rank[k]++;//֧���k����ĸ���+1
				if (offspring[k] < offspring[j] && !(offspring[k] == offspring[j]))//k֧��j
				{
					offspring[k].sp++;//k��֧��ļ���+1
					int m = offspring[k].sp - 1;
					cset[k][m] = j;//cset[2N][2N],cset[k][m]=j��ʾ��k������֧��ĵ�m������j
				}
			}
		}
		printf("֧���%d����ĸ���Ϊ%d��\n", k, rank[k]);
		printf("��%d��֧��ĵ�Ϊ��", k);
		cout << offspring[k].sp << "\n";
	}
	int curr_rank = 0;
	while (1)
	{
		int stop_count = 0;
		int* rank2 = new int[size];
		for (int k = 0; k < size; k++)
			rank2[k] = rank[k];//����
		for (int k = 0; k < size; k++)
		{
			if (offspring[k].np == -1 && rank[k] == 0)
			{
				offspring[k].np = curr_rank; //��֧�����Ϊ0�ĵ��������0Ϊ���
				for (int j = 0; j < offspring[k].sp; j++)//��k��֧��ļ�����
				{
					int id = cset[k][j];//����k��֧��ĵ�j+1�����id
					rank2[id]--; //��֧��id�ĵ�ø�����һ
					stop_count++;
				}
			}
		}
		for (int k = 0; k < size; k++)
		{
			rank[k] = rank2[k];//�����º�֧��k�ĵ������ƻ�rank[]
		}
		delete[] rank2;
		curr_rank++;
		//cout<<"curr_rank: "<<curr_rank<<"\n\n\n";
		if (stop_count == 0)
			break;
	}
	delete[] rank;
	for (int i = 0; i < size; i++)
		delete cset[i];
	delete[] cset;
}

//ӵ��������ѡ�񣬴˹����ڷ�֧������֮��
void zc_NSGA2::evaluation_density()
{
	population.clear();
	int size = offspring.size();     //size = 2N
	int rank = 0;
	while (1)
	{
		int count = 0;
		for (int i = 0; i < size; i++)
		{
			if (offspring[i].np == rank)
				count++;
		}
		//cout<<"count: "<<count<<"\n\n";
		int size2 = population.size() + count; //���ŵ�����
		//cout<<"���ŵ�����size2: "<<size2<<"\n\n";
		if (size2 > pops)
			break;
		for (int i = 0; i < size; i++)
			if (offspring[i].np == rank)
				population.push_back(offspring[i]);//�������и��壬�ҵ�����֧��ĵ�
		rank++;
		//cout<<"rank:"<<rank<<"\n\n";
		if (population.size() >= pops)
			break;
	}//����֧��ⰴ����̶Ȳ�����䵽population�У�ֱ����С����pops����Ϊֹ,�����ʱ����������ӵ���Ⱦ���������
	if (population.size() < pops)
	{
		vector <zc_NSGA2Ind> list;//��������������һ���ĸ���
		for (int i = 0; i < size; i++)
			if (offspring[i].np == rank)
				list.push_back(offspring[i]);
		int s2 = list.size(); //������pops���Ǹ������ĵ������
		//cout<<"����pops���Ǹ������ĵ������s2: "<<s2<<"\n\n";
		double *density = new double[s2];
		int *idx = new int[s2];
		int *idd = new int[s2];
		double *obj = new double[s2];
		for (int i = 0; i < s2; i++)
		{
			density[i] = 0;
			idx[i] = i;
		}//��ʼ��
		for (int j = 0; j < nobj; j++)
		{
			for (int i = 0; i < s2; i++)
			{
				idd[i] = i;
				obj[i] = list[i].y_obj[j];
			}
			minfastsort(obj, idd, s2, s2);
			density[idd[0]] += -1.0e+30;//�߽�ֵȡ�����
			density[idd[s2 - 1]] += -1.0e+30;
			for (int k = 1; k < s2 - 1; k++)
			{
				density[idd[k]] += -(obj[k] - obj[k - 1] + obj[k + 1] - obj[k]);//��������ӵ���Ⱦ���								
				//cout<<"density[] = "<<density[idd[k]]<<"\n\n";
			}
		}
		delete[] idd;
		delete[] obj;
		int s3 = pops - population.size();
		//	cout<<"s3 = "<<s3<<"\n\n";
		minfastsort(density, idx, s2, s3);
		for (int i = 0; i < s3; i++)
			population.push_back(list[idx[i]]);
		delete[] density;
		delete[] idx;
	}
	offspring.clear();
}

//��N������ŵ�offspring�У�����֤û���ظ�
void zc_NSGA2::fill_union(zc_NSGA2Ind &ind)
{
	bool flag = true;
	int size = offspring.size();//����offspring�и��������
	//cout<<"����offspring�и��������: "<<size<<"\n\n";
	for (int i = 0; i < size; i++)
		if (ind == offspring[i])
		{
			flag = false;
			break;
		}
	if (flag)
		offspring.push_back(ind);
}

//�����鰴�մ�С��������
void minfastsort(double x[], int idx[], int n, int m)
{
	for (int i = 0; i < m; i++)
	{
		for (int j = i + 1; j < n; j++)
			if (x[i] > x[j])
			{
				double temp = x[i];
				x[i] = x[j];
				x[j] = temp;
				int temp2 = idx[i];
				idx[i] = idx[j];
				idx[j] = temp2;
			}
	}
}

/*static unsigned int my_int = 1;
double rnd_uni(long *l)
{
	//srand((unsigned int)*l);
	//srand(time(NULL));
	my_int = (my_int + 500) % (INT_MAX - 50);
	srand(my_int);
	return rand() * 1.0 / RAND_MAX;
}*/

union FloatRand
{
	struct
	{
		unsigned long Frac : 23;
		unsigned long Exp : 8;
		unsigned long Signed : 1;
	} BitArea;
	float Value;
	unsigned long Binary; /* for debug only */
};

double rnd_uni(long *l)
{
	/*union FloatRand r;

	r.BitArea.Signed = 0;
	r.BitArea.Exp = 1;
	r.BitArea.Frac = (rand() * rand()) % 0x800000;
	if (r.BitArea.Frac == 0x7FFFFF)
		r.BitArea.Exp = 0x7D;
	else if (r.BitArea.Frac == 0)
		r.BitArea.Exp = 0x7E;
	else
		r.BitArea.Exp = 0x7E;

	//return r.Value;
	return (double)(r.Value - 0.5)*2.0;*/

	return rand() * 1.0 / RAND_MAX;
}