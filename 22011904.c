#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define MAX 100

void fnc_call(double polynomial[2][MAX],int count)
{

    int i;
    for(i=0; i<count; i++){
        printf("\nPolynomial: x_coef * x ^ x_exp");
        printf("\nx_coef: ");
        scanf("%lf", &polynomial[0][i]);
        printf("x_exp: ");
        scanf("%lf", &polynomial[1][i]);
        printf("Added: %lf * x^%lf", polynomial[0][i], polynomial[1][i]);
    }
    printf("\nPolynom is: ");
    for(i=0; i<count; i++){
        printf(" (%lf * x^%lf) + ",polynomial[0][i], polynomial[1][i]);
    }

}


double f(double polynomial[2][MAX], double x, int count)
{
    int i;
    double result = 0;
    for(i=0; i<count; i++){
        result += (polynomial[0][i]*pow(x,polynomial[1][i]));
    }
    return result;
}

double df(double polynomial[2][MAX], int count,  double poly_drv[2][MAX], double x)
{
    int	i;
    double fnc_derivative = 0;
	for(i=0; i<count; i++){
		poly_drv[0][i] = polynomial[0][i];
		poly_drv[1][i] = polynomial[1][i];
	}

	for(i=0; i<count; i++){
		if(poly_drv[1][i] == 0){
			poly_drv[1][i] = 0;
		}else{
			poly_drv[0][i] *= poly_drv[1][i];
			poly_drv[1][i]--;
		}
	}
	for(i=0; i<count; i++){
        fnc_derivative += (poly_drv[0][i]*pow(x,poly_drv[1][i]));
    }
    return fnc_derivative;
}

void bisection_method()
{
    double	polynom[2][MAX] = {0};
    double a, b, c,c_old, epsilon,error;
    int iter = 1, MAX_ITER, count;

    printf("Please enter polynomial count: ");
    scanf("%d",&count);
    fnc_call(polynom, count);
    printf("\nEnter the interval [a,b]: ");
    scanf("%lf %lf", &a, &b);

    if(f(polynom,a,count)*f(polynom,b,count) < 0){
    printf("Enter MAX number of iterations and epsilon: ");
    scanf("%d %lf", &MAX_ITER, &epsilon);

    do{
        c_old = c;
        c = (a + b) / 2.0;
        error = abs(c-c_old);
        printf("\na = %lf\nb = %lf\nc = %lf\niteration = %d\nea = %lf\n", a,b,c,iter,error);
        if(error < epsilon){
            printf("\nSolution found: x = %lf\n", c);
            main();
        }

        if(f(polynom,a,count)*f(polynom,c,count) < 0){b = c;}
        else{a = c;}
        iter++;
    }while(((b-a)/pow(2,iter)) > epsilon);
    if(iter>MAX_ITER)printf("Error: Maximum iterations exceeded\n");
    main();
    }
    else printf("f(a) and f(b) have the same sign\n");
    main();
}

void regula_falsi()
{
    double a, b, c, fa, fb, fc, epsilon;
    double polynom[2][MAX];
    int i = 0, iter, count;
    printf("Please enter polynomial count: ");
    scanf("%d",&count);
    fnc_call(polynom, count);
    printf("\nEnter the initial interval [a, b]: ");
    scanf("%lf %lf", &a, &b);
    printf("Enter the maximum number of iterations: ");
    scanf("%d", &iter);
    printf("Enter the error tolerance: ");
    scanf("%lf", &epsilon);
    fa = f(polynom,a,count);
    fb = f(polynom,b,count);
    while(i < iter || ((b-a)/pow(2,iter)) < epsilon) {
        c = (a*fb - b*fa)/(fb - fa);
        fc = f(polynom,c,count);
        printf("Iteration %d\na = %lf\nf(a) = %lf\nb = %lf\nf(b) = %lf\nc = %lf\nf(c) = %lf\n", i+1,a,fa,b,fb,c,fc);
        if(abs(fc) < epsilon){
            printf("Solution found: x = %lf\n", c);
            main();
            }
        else if((fa*fc) < 0){
                b = c;
                fb = fc;}
        else{
            a = c;
            fa = fc;}
            i++;}
    printf("Error: Regula-Falsi algorithm did not converge\n");
    main();
}

void newton_raphson()
{
    double x0, x1, epsilon, fx0, df_x0,polynom[2][MAX],poly_drv[2][MAX];
    int i = 0, iter,count;
    printf("Please enter polynomial count: ");
    scanf("%d",&count);
    fnc_call(polynom, count);
    printf("\nEnter the initial guess: ");
    scanf("%lf", &x0);
    printf("Enter the maximum number of iterations: ");
    scanf("%d", &iter);
    printf("Enter the error tolerance: ");
    scanf("%lf", &epsilon);
    while (i < iter){
        fx0 = f(polynom,x0,count);
        df_x0 = df(polynom,count-1,poly_drv,x0);
        if(df_x0 == 0){
            printf("Derivative of function equals to 0...");
            main();
        }
        x1 = x0 - (fx0/df_x0);
        printf("Iteration %d\nx = %lf\nf(x) = %lf\n", i+1, x1, f(polynom,x1,count));

        if (abs(x1 - x0) < epsilon) {
            printf("Solution found: x = %lf\n", x1);
            i = iter;
             main();
        }
        x0 = x1;
        i++;
    }
    printf("Error: Newton-Raphson algorithm did not converge\n");
    main();
}



void inverse()
{
    int i,j,k,N;
    double A[MAX][MAX], B[MAX][MAX], diagonal, factor;
    printf("Enter max number of rows and columns: ");
    scanf("%d", &N);
    printf("Enter elements of matrix A: ");
    for(i=0; i<N; i++){
        for(j=0; j<N; j++){
            scanf("%lf", &A[i][j]);
        }
    }
    printf("Matrix A:\n");
    for(i=0; i<N; i++){
        for(j=0; j<N; j++){
            printf("%lf ", A[i][j]);
        }printf("\n");
    }
    for(i=0; i<N; i++){
        for(j=0; j<N; j++){
            if(i==j){B[i][j] = 1;}
            else B[i][j] = 0;
        }
    }
    for(i=0; i<N; i++){
		diagonal = A[i][i];
		for(j=0; j<N; j++){
			A[i][j] /= diagonal;
			B[i][j] /= diagonal;
		}
		for(j=0; j<N; j++){
			if(i!=j){
				factor = A[j][i];
				for(k=0; k<N; k++){
					A[j][k] = A[j][k] - A[i][k]*factor;
					B[j][k] = B[j][k] - B[i][k]*factor;
				}
			}
		}
	}
	printf("\nInverse matrix is:\n");
	for(i=0; i<N; i++){
		for(j=0; j<N; j++){
			printf("%lf ",B[i][j]);
		}printf("\n");
	}
    main();
}
void gauss_elimination()
{
    double A[MAX][MAX+1];
    double x[MAX],f;
    int N,i,j,k,l;
    printf("Enter size of matrix: ");
    scanf("%d", &N);
    printf("Enter the matrix A|B:\n");
    for (i=0; i<N; i++) {
        for (j=0; j<N+1; j++) {
            scanf("%lf", &A[i][j]);
        }
    }
    printf("Matrix before the elimination: \n");
    for(i=0; i<N; i++){
        for(j=0; j<N; j++){
            printf("%d  ", A[i][j]);
        }printf("\n");
    }
    for(i=0; i<N-1; i++) {
        for(j=i+1; j<N; j++) {
            f = A[j][i] / A[i][i];
            for(k=i; k<N+1; k++) {
                A[j][k] = A[j][k] - f * A[i][k];
            }
        }
        printf("Matrix after the %dth step of forward elimination:\n", i+1);
        for(k=0; k<N; k++) {
            for(l=0; l<N+1; l++) {
                printf("%lf ", A[k][l]);
            }
            printf("\n");
        }
        printf("\n");
    }
    for(i=N-1; i>=0; i--) {
        x[i] = A[i][N];
        for(j=i+1; j<N; j++) {
            x[i] = x[i] - A[i][j] * x[j];
        }
        x[i] = x[i] / A[i][i];
    }
    printf("Resultant matrix:\n");
    for(i=0; i<N; i++) {
        for(j=0; j<N+1; j++) {
            printf("%lf ", A[i][j]);
        }
        printf("\n");
    }
   main();
}


void gauss_seidal()
{
    int i,j,n, max_iterations, iter = 0;
    double A[MAX][MAX], b[MAX], x[MAX], old_x[MAX], ea = 0.0, epsilon, tmp;
    printf("Enter the size of the system: ");
    scanf("%d", &n);
    printf("Enter the elements of the matrix:\n");
    for(i=0; i<n; i++){
        for(j=0; j<n; j++) {
            scanf("%lf", &A[i][j]);
        }
    }
    printf("Enter the elements of the vector (x1,x2,x3...):\n");
    for(i=0; i<n; i++){
        printf("x%d = ", i+1);
        scanf("%lf", &b[i]);
    }
        for(i=0; i<n; i++) {
        x[i] = 0.0;
    }
    printf("Enter the maximum number of iterations: ");
    scanf("%d", &max_iterations);
    printf("Enter the tolerance: ");
    scanf("%lf", &epsilon);
    iter = 0;
    do{
        for(i=0; i<n; i++){
            old_x[i] = x[i];
            tmp = 0.0;
            for(j=0; j<n; j++) {
                 if(i!=j){
                    tmp+=(A[i][j]*x[j]);
                }
            }
            x[i] = (b[i] - tmp) / A[i][i];
            ea += abs(x[i] - old_x[i]);
            if(ea < epsilon){
                iter = max_iterations;
            }
        }
        iter++;
    } while(iter < max_iterations);
    if(ea < epsilon){
        printf("Solution vector: \n");
        for(i=0; i<n; i++) {
            printf("x%d = %lf\n", i+1, x[i]);
        }
        main();
    } else{
        printf("Gauss-Seidel iteration failed to converge after %d iterations.\n", max_iterations);
    }
    main();
}


void simpson()
{
    double a,b,h,sum,x,polynom[2][MAX];
    int n,i,count;
    double x0,x1,x2,fx0,fx1,fx2, I = 0.0;
    printf("Please enter polynomial count: ");
    scanf("%d",&count);
    fnc_call(polynom, count);
    printf("\nEnter the interval [a;b] : ");
    scanf("%lf%lf", &a, &b);
    printf("\nEnter the number of sub-intervals: ");
    scanf("%d", &n);
    if((n % 2) != 0){
        printf("\nError! Number of sub-intervals must be EVEN!\nEnter again: ");
        scanf("%d" ,&n);
    }
    h = (b - a)/n;
    x0 = a;
    x1 = a + h;
    x2 = a + 2*h;
    for(i=0; i<n/2; i++) {
        fx0 = f(polynom,x0,count);
        fx1 = f(polynom,x1,count);
        fx2 = f(polynom,x2,count);
        I += (h/3) * (fx0 + 4*fx1 + fx2);
        x0 = x2;
        x1 = x2 + h;
        x2 = x2 + 2*h;
    }
    printf("\nArea of function f(x) = %lf\n", I);
    main();
}


void derivative()
{
    double x, h, cent_derv,polynom[2][MAX];
    int count;
    printf("Please enter polynomial count: ");
    scanf("%d",&count);
    fnc_call(polynom, count);
    printf("\nEnter the value of x: ");
    scanf("%lf", &x);
    printf("Enter the value of h: ");
    scanf("%lf", &h);
    cent_derv =  (f(polynom,x + h,count) - f(polynom,x - h,count)) / (2 * h);
    printf("\nDerivative of %lf at %lf is %lf\n", x,h,cent_derv);
    main();

}


void trapezoid()
{
    double a,b,h,n,sum,x,area,polynom[2][MAX];
    int i,count;
    printf("Please enter polynomial count: ");
    scanf("%d",&count);
    fnc_call(polynom, count);
    printf("\nEnter the interval [a;b]: ");
    scanf("%lf%lf", &a, &b);
    printf("\nEnter the number of sub-intervals: ");
    scanf("%lf", &n);
    h = (b-a)/n;
    x = a;
    sum = (f(polynom,a,count)+f(polynom,b,count))/2;

    for(i=1; i<n; i++) {
        x += h;
        sum += f(polynom,x,count);
    }
    area = h*sum;
    printf("\nArea of function f(x) = %lf\n", area);
        main();

}


void gregory_newton()
{
    int i,j,n,count;
    double x[MAX], polynom[2][MAX], difference[MAX][MAX];
    double tmp, result, value;
    printf("Please enter polynomial count: ");
    scanf("%d", &count);
    fnc_call(polynom, count);

    printf("\nEnter the number of data points: ");
    scanf("%d", &n);
    for(i=0; i<n; i++) {
        printf("Enter x[%d]: ",i);
        scanf("%lf", &x[i]);
    }
    printf("Enter the value of x at which to evaluate the polynomial: ");
    scanf("%lf", &value);

    result = 0;
    for(i=0; i<n; i++) {
        difference[i][0] = f(polynom,x[i],count);
    }
    for(j=1; j<n; j++) {
        for(i=j; i<n; i++) {
            difference[i][j] = (difference[i][j-1] - difference[i-1][j-1]) / (x[i] - x[i-j]);
        }
    }
    for(j=1; j<n; j++) {
        tmp = difference[j][j];
        for(i=0; i<j; i++) {
            tmp *= (value - x[i]);
        }
        result += tmp;
    }
    printf("The value of the polynomial at x = %lf is: %lf\n", value, result);
    main();
}

int main()
{
    int choice;
    printf("Choose from the methods below(1-10);\n1. Bisection method\n2. Regula-Falsi method\n3. Newton-Rapshon method\n4. Inverse of a matrix(N*N)\n5. Gauss Elimination\n6. Gauss-Seidal method\n7. Numerical Derivative\n8. Simpson Method\n9. Trapezoid method\n10. Gregory-Newton\nPress 0 to EXIT.");
    scanf("%d", &choice);
    do{
    switch(choice)
    {
    case 1:
        bisection_method();
        break;
   case 2:
        regula_falsi();
        break;
    case 3:
        newton_raphson();
        break;
    case 4:
        inverse();
        break;
    case 5:
        gauss_elimination();
        break;
    case 6:
        gauss_seidal();
        break;
    case 7:
        derivative();
        break;
    case 8:
        simpson();
        break;
    case 9:
        trapezoid();
        break;
    case 10:
        gregory_newton();
        break;
    default:
        break;
    }

    }while(choice != 0);

    return 0;
}
