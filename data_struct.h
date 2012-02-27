

//==============================================================================================//
//			DEFINICION DE LAS ESTRUCTURAS DE DATOS					//
//==============================================================================================//


struct genSenal	{ 
		float afinidad[MAXGEN_REG][3];
		float produccion;
		float coef_difusion;
		float difusion;
		};

struct genRegulador
		{ 
		float afinidad[MAXGEN_SEN + MAXGEN_EST][3];
		float degradacion;
		float produccion;
		float actividad;
		};

struct genEstructural
		{ 
		//float afinidad[MAXGEN_SEN + MAXGEN_FUN];
		char tipo;
		float actividad;
		float produccion;
		float degradacion;
		};

struct Coord	{
		int x;
		int y;
		int z;
		};

struct propFis	{
		float k;
		float mu;
		float l0;
		int vecinos;
		};


//==============================================================================================//
//			ESTRUCTURA DE DATOS PRIMARIA: celula					//
//==============================================================================================//


struct celula	{
		struct genSenal 	Senales[MAXGEN_SEN];
		struct genRegulador 	Reguladores[MAXGEN_REG];
		struct genEstructural 	Estructurales[MAXGEN_EST];
		struct Coord 		posicion;
		struct propFis 		fisica;
		int vida_muerte;
		};

struct celula *crear_celula(int x, int y, int z);
struct celula *celulas[MAX_CEL];
int celulas_TOT=0;
int inicializado=0;
int aff_senales[MAXGEN_SEN][MAXGEN_REG][3];
int aff_reguladores[MAXGEN_REG][MAXGEN_SEN+MAXGEN_EST][3];

//==============================================================================================//
//			OTROS DATOS GLOBALES NECESARIOS PARA EL PROGRAMA			//
//==============================================================================================//

int espDifusivo2D[LARGO][ANCHO][MAXGEN_SEN+1];
int espDifusivo2D_old[LARGO][ANCHO][MAXGEN_SEN];
float espFisico[1][1][1];
float t_simulacion;
