/* --------------------------------------------------*
 |	 Código de ferramentas para trabalhar com cônicas|
 | 	 Feito por:                                      |
 |		Lourenço de Salles Roselino - 11796805       |
 |		Milena Corrêa da Silva - 11795401            |
 |	Para a materia de Geometria Analitica SMA0300    |
 |	BCC 020 - ICMC                                   |
 * --------------------------------------------------*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define AMOUNT_OF_COEFFICIENTS 6
#define ROUND_NUMBER 10
#define TRUE 1
#define FALSE 0

typedef enum _coefficients {
    A, // quadrático (x2)
    B, // quadratico misto (xy)
    C, // quadrático (y2)
    D, // linear (x)
    E, // linear (y)
    F, // termo independente
} Coefficients;

/* Retorna o minimo multiplo comum de dois número (necessario para exibir na forma fracionada) */
int gcd(int n1, int n2) {
	int remainder;
	while (n2) {
		remainder = n1 % n2;
		n1 = n2;
		n2 = remainder;
	}
	return n1;
}

/* Verifica se o double tem numeros depois da virgula */
int isRational(double n) {
	if(n - floor(n) != 0) {
		return TRUE;
	}
	return FALSE;
}


/* Elimina os termos lineares a partir da translação */
void translation(double *coefficients) {
    double x,y;
    /* Calcula os determinantes que seram usados na regra de Cramer a partit da 1 e 2 linhas da matriz simétrica g */
    double mainDeterminant = (coefficients[A] * coefficients[C]) - (pow(coefficients[B]/2, 2));
    double xDeterminant = -(coefficients[D]/2 * coefficients[C]) + (coefficients[B]/2 * coefficients[E]/2);
    double yDeterminant = -(coefficients[E]/2 * coefficients[A]) + (coefficients[B]/2 * coefficients[D]/2);
    

    if (mainDeterminant == 0 && ((xDeterminant != 0) | (yDeterminant != 0)))  { 
        /* Caso o sistema for impossível */
        printf("Function impossible to translate\n");
        return;
    }

    if (mainDeterminant == 0 && xDeterminant == 0 && yDeterminant == 0) { 
        /* Caso o sistema for possível mas indeterminado. Assim, assumimos um número x qualquer e calculamos o y a partir desse x */
        x = 0; 
        y = -(coefficients[D]/coefficients[B]);
    }

    else { // Sistema possível e determinado. Calculamos o x e y a partir da regra de Cramer
        x = xDeterminant / mainDeterminant;
        y = yDeterminant / mainDeterminant;
    }

    /* Usamos a 3 linha da matriz simetrica g */
    double translatedIndependentTerm = (coefficients[D]/2 * x) + (coefficients[E]/2 * y) + coefficients[F];

    /* Mudamos os coefiientes */ 
    coefficients[D] = coefficients[E] = 0;
    coefficients[F] = translatedIndependentTerm;
}

/* Elimina o termo quadŕatico misto a partir da rotação */ 
void rotation(double *coefficients) {
    double cotan2 = (coefficients[A] - coefficients[C]) / coefficients[B];
    double sine2 = 1 / sqrt(1 + pow((coefficients[A]-coefficients[C])/coefficients[B], 2));
    double cosine2 = cotan2 * sine2;
    /* Novamente, calcula os determinantes que seram usados na regra de Cramer */
    double mainDeterminant = -2.0; // Como o determinante principal nunca sera 0, o sistema é sempre determinado
    double _aDeterminant = -(coefficients[A] + coefficients[C]) -(coefficients[B]/sine2);
    double _cDeterminant = -(coefficients[A] + coefficients[C]) + (coefficients[B]/sine2);
    

    double sine = sqrt((1 - cosine2) / 2);
    double cosine = sqrt((cosine2 + 1) / 2 ); 
    double tmpD = coefficients[D];
    double tmpE = coefficients[E];

    if (coefficients[D]) 
        coefficients[D] = (tmpD * cosine) + (tmpE * sine);
    if (coefficients[E])
        coefficients[E] = (cosine * tmpE) - (sine * tmpD);
    

    /* Efetiva a regra de Cramer e muda os  coeficientes */
    coefficients[A] = _aDeterminant / mainDeterminant;
    coefficients[C] = _cDeterminant / mainDeterminant;
    coefficients[B] = 0;

}

/* Identifica o tipo de cônica a partir invariantes ortogonais da conica*/
char *getConicType(double *coeficients) {
    double discriminant1 = coeficients[A] + coeficients[C]; // traço de M2 excluindo b
    double discriminant2 = coeficients[A] * coeficients[C]; // determinante de M2 excluindo b
    double discriminant3 = (coeficients[A] * coeficients[C] * coeficients[F]) -  
    ((pow(coeficients[D], 2 ) * coeficients[C] + pow(coeficients[E], 2 ) * coeficients[A]) / 4) ; //determinante de M3 excluindo b
    double delta = (coeficients[A] * coeficients[F]) + (coeficients[F] * coeficients[C]); //cofatores

    if (discriminant2 != 0) {
        if (discriminant3 != 0) {
            if (discriminant2 < 0) 
                return "hiperbole";
            if (discriminant1 * discriminant3 > 0) 
                return "conjunto vazio";

            if (coeficients[A] != coeficients[C] || coeficients[B] != 0) 
                return "elipse";

            return "circunferencia";
    
            }        

        if (discriminant2 < 0) 
            return "duas retas concorrentes";

        return "ponto";
    }

    if(discriminant3 != 0) 
        return "parabola";
    if (delta < 0)
        return "duas retas paralelas";
    if(delta == 0) 
        return "reta";
    
    return "conjunto vazio";
    
}

/* Vai criar um vetor de strings com as variaveis na ordem a/fx² + b/fxy + c/fy² + dx + ey + f = 0 */
char** getResultString(double* coefficients) {
	char** result;
	result = malloc(sizeof(char*) * AMOUNT_OF_COEFFICIENTS);
	int buffer;
	for(int i = 0; i < AMOUNT_OF_COEFFICIENTS; i++) {
		result[i] = malloc(sizeof(char) * 4096);
		if(isRational(coefficients[i])) {
			buffer = round(coefficients[i] * ROUND_NUMBER);
			sprintf(result[i], "(%d/%d)", buffer/gcd(buffer, ROUND_NUMBER), ROUND_NUMBER/gcd(buffer, ROUND_NUMBER));
		} else {
			sprintf(result[i], "%d", (int) coefficients[i]);
		}
		result[i] = realloc(result[i], strlen(result[i]));
	}

	return result;
}


int main() {
    double tmpNum;
	double angle;
    double *coefficients = malloc(sizeof(double) * AMOUNT_OF_COEFFICIENTS);
	char** results;
    for (int i = 0; i < AMOUNT_OF_COEFFICIENTS; ++i) {
        scanf("%lf", &tmpNum);
        coefficients[i] = tmpNum;
    }


    translation(coefficients);
    for (int i = 0; i < AMOUNT_OF_COEFFICIENTS; ++i) {
        printf("%.1lf ", coefficients[i]);
    }
    putchar('\n');


    if (coefficients[B]) {
		angle = atan(coefficients[B] / coefficients[A] * coefficients[C] ) / 2;
		printf("%lf\n", coefficients[B] / coefficients[A] * coefficients[C]);
		printf("θ = %lf rads\n", angle);
		rotation(coefficients);
	}
	results = getResultString(coefficients);

	if(coefficients[A]) {
		printf("%sx²", results[A]);
	}
	if(coefficients[B]) {
		printf(" + %sxy", results[B]);
	}
	if(coefficients[C]) {
		printf(" + %sy²", results[C]);
	}
	if(coefficients[D]) {
		printf(" + %sx", results[D]);
	}
	if(coefficients[E]) {
		printf(" + %sy", results[E]);
	}
	if(coefficients[F]) {
		printf(" + %s", results[F]);
	}
	printf(" = 0\n");

    printf("Tipo de conica: %s\n", getConicType(coefficients));

	if(!strcmp(getConicType(coefficients),"elipse")) {
		double p = sqrt(fabs((coefficients[F] / coefficients[A]) - (coefficients[F] / coefficients[C]))); // Calculo do valor dos focos, sendo p o valor focal
		printf("p = %lf\n", p);
		printf("Focos: F1=(%lf, %lf), F2=(%lf, %lf)\n", coefficients[D] + p, coefficients[E], coefficients[D] - p, coefficients[E]);
		printf("Centro: C=(%lf,%lf)\n", coefficients[D]/coefficients[F], coefficients[E]/coefficients[F]);
		printf("Vertices: V1 = (±√%lf, 0) V2 = (0, ±√%lf)\n", coefficients[F]/coefficients[A], coefficients[F]/coefficients[C]); // Isso vai dar um resultado incorreto se não tivermos eliminado os termos lineares, depois eu resolvo isso melhor
	}

	for(int i = 0; i<AMOUNT_OF_COEFFICIENTS; i++) {
		free(results[i]);
	}
	free(results);
    free(coefficients);

}
