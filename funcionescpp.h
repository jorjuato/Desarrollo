

int configura(void);

int inicializa(void);

int simula(void);

int difusion_engine(void);

int fisica_engine(void);

int red_genetica_engine(void);

int genes_estructurales(int i);

int genes_reguladores(int i);

int genes_senales(int i);

int genes_fisica(int i);

int muerte(int i);

struct celula *crear_celula(int x,int y,int z);

int search(int *coor, int i);

int retornar_indice(int x, int y, int z);

int visualiza(void); 

int print_init_conf(void);

int print_red_genetica(void);

int print_espDifusivo(int);

int print_genestructurales(void);

int vidaymuerte(void);

int dibujar_red_genetica(void);

int dibujar_difusion(void);









