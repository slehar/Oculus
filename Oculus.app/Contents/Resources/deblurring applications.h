typedef struct{
	int inputint1;
	int inputint2;
	int inputint3;
	int inputint4;
	float inputfloat1;
	float inputfloat2;
	float inputfloat3;
	float inputfloat4;
	float inputfloat5;
	float inputfloat6;
	float inputfloat7;
	float inputfloat8;
	float inputfloat9;
	float inputfloat10;
	float inputfloat11;
	float inputfloat12;
	float inputfloat13;
	float inputfloat14;
}DEBLURINPUT;

void copyToWindowArray(unsigned char *,unsigned char *,unsigned char *,int, char *,int,int,int);
void processImage(short,short,short,char *,char *,DEBLURINPUT *);