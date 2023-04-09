#include <stdio.h> 
#include <math.h>

double F(double x0, double x1)
{
    return 100*(x1-x0*x0)*(x1-x0*x0) + (1-x0)*(1-x0);
}

int write_file(char *filename)
{
    int i;
    double x0, x1, x, delta;
    x0 = -2;
    x1 = 1;
    delta = 0.04;
    FILE *p_file; 
    
    p_file = fopen(filename, "w");
    for (i = 0; i <= 100; i++)
    {
        x = x0 + i*delta;
        fprintf(p_file, "%e, %e\n", x, F(x, x1));
    }
    fclose(p_file);
    return 0;
}

typedef struct pcoordinate //initialise p coordinates
{
    double position[2];
} pcoordinate;

double pbarf(double vec1[2], double vec2[2], int N) //adds two vectors together then divides by 2 (finds midpoint/ centroid), 'f' at end indicates that it's a function
{
    double result[2];
    int i;
    for (i = 0; i < 2; i++)
    {
        result[i] = (vec1[i] + vec2[i]) / 2;
    }

    return result[N]; // returns nth element of pbar vector. Note that N = 0 corresponds to first element
    
}

double ybarf(double y1, double y2, double y3)
{
    double result;
    result = y1+y2+y3;
    return result;
}

double pstarf(double vec1[2], double vec2[2], int N) //calculates p* transformtion (multiplies pbar by 2 then subtracts p_l)
{
    double result[2];
    int i;
    for (i = 0; i < 2; i++)
    {
        result[i] = 2*vec1[i] - vec2[i];
    }
    
    return result[N];
} 

double pstarstar1f(double vec1[2], double vec2[2], int N)
{
    double result[2];
    int i;
    for (i = 0; i < 2; i++)
    {
        result[i] = 2*vec1[i] - vec2[i];
    }
    return result[N];
}

double pstarstar2f(double vec1[2], double vec2[2], int N)
{
    double result[2];
    int i;
    for (i = 0; i < 2; i++)
    {
        result[i] = (vec1[i] + vec2[i]) / 2;
    }
    return result[N];
}

double rosenbrockf(double vector[2]) //function finds 'elevation' of coordinate
{

    double y;
    y = 100*(vector[1] - vector[0]*vector[0])*(vector[1] - vector[0]*vector[0]) + (1-vector[0])*(1-vector[0]);
    return y;
} 

double standevf(double y1, double y2, double y3, double ybar) // function finds if 'terminating' condition is satisfied
{
    double result = 0;
    if ( sqrt( (pow(y1-ybar, 2) / 2) + (pow(y2 - ybar, 2) / 2) + (pow(y3 - ybar,2) / 2) ) > pow(10, -8))
    {
        result = 1;
    }
    return result;
}


int algorithm() //contains conditional statements that result in answer: bulk of program
{
    pcoordinate p1, p2, p3, pl, ph, pm, pbar, pstar, pstarstar; //sets initial coordinates and defines variables
    p1.position[0] = 0.0;
    p1.position[1] = 0.0;

    p2.position[0] = 2.0;
    p2.position[1] = 0.0;

    p3.position[0] = 0.0;
    p3.position[1] = 2.0;


    int i, i2; //i is used in all for loops, i2 is used in nested for loops
    
    double y1, y2, y3, ybar; // uses y1 (the number) not to be confused with yl (the letter)!

    double count;
    double *pcount;
    pcount = &count; //used a pointer to count the number of iterations the program undergoes

    
    double y_l, y_m, y_h, ystar, ystarstar; //uses y_l the letter!

    //looping algorithm starts
    do
    {

    *pcount = *pcount + 1;
    
    // calculate elevation of each coordinate
    y1 = rosenbrockf(p1.position);
    y2 = rosenbrockf(p2.position);
    y3 = rosenbrockf(p3.position);

    //better sorting


    
    //sorts each y elevation, relates this to it's corresponding p coordinate
    if (y1 > y2)
    {
        if (y3 > y2)
        {
            if (y3 > y1)
            {
                for (i = 0; i < 2; i++)
                {
                    ph.position[i] = p3.position[i];
                    pm.position[i] = p1.position[i];
                    pl.position[i] = p2.position[i];
                }
                // y3 biggest, y1 in middle, y2 lowest
            }
            else
            {
                for (i = 0; i < 2; i++)
                {
                    ph.position[i] = p1.position[i];
                    pm.position[i] = p3.position[i];
                    pl.position[i] = p2.position[i];
                }
                //y1 biggest, y3 mid, y2 low
            }
        }
        else
        {
            for (i = 0; i < 2; i++)
                {
                    ph.position[i] = p1.position[i];
                    pm.position[i] = p2.position[i];
                    pl.position[i] = p3.position[i];
                }
            //y1 big, y2 mid, y3 low
        }
        
    }
    else
    {
        if (y3 > y2)
        {
            for (i = 0; i < 2; i++)
                {
                    ph.position[i] = p3.position[i];
                    pm.position[i] = p2.position[i];
                    pl.position[i] = p1.position[i];
                }
            // y3 big, y2 mid, y1 smol
        }
        else
        {
            if (y1 > y3)
            {
                for (i = 0; i < 2; i++)
                {
                    ph.position[i] = p2.position[i];
                    pm.position[i] = p1.position[i];
                    pl.position[i] = p3.position[i];
                }
                // y2 big, y1 mid, y3 smol
            }
            else
            {
                for (i = 0; i < 2; i++)
                {
                    ph.position[i] = p2.position[i];
                    pm.position[i] = p3.position[i];
                    pl.position[i] = p1.position[i];
                }
                // y2, y3, y1
            }
        } 
        
    }

    y_l = rosenbrockf(pl.position);
    y_m = rosenbrockf(pm.position);
    y_h = rosenbrockf(ph.position);

    // sorting complete!
    // Find pbar, pstar, ystar
    for (i = 0; i < 2; i++) //finds pbar
    {
        pbar.position[i] = pbarf(pl.position, pm.position, i);
    }

    for (i = 0; i < 2; i++)
    {
        pstar.position[i] = pstarf(pbar.position, ph.position, i);
    }

    ystar = rosenbrockf(pstar.position);
    
    int a,b,c,d,e,f; // these will represent each 'decision' in the flow chart
    a=0;
    b=0;
    c=0;
    d=0;
    e=0;
    f=0;
    
    // is y* < y_l? Decision 'a'
    if (ystar < y_l) 
    {
        a = 1;
    }
    
    if (a==1) //Decision 'd'
    {
        for (i = 0; i < 2; i++) //calculate pstarstar
        {
            pstarstar.position[i] = pstarstar1f(pstar.position, pbar.position, i);
        }
        
        ystarstar = rosenbrockf(pstarstar.position); //calculates ystarstar
        
        if (ystarstar < y_l)
        {
            for (i = 0; i < 2; i++) //replace ph by pstarstar
            {
                ph.position[i] = pstarstar.position[i];
            }
            
        }
        else
        {
            for (i = 0; i < 2; i++) //replace ph by pstar
            {
                ph.position[i] = pstar.position[i];
            }
            
        }
    }

    if (a==0) //Decision 'b'
    {
        if (ystar > y_m)
        {
            b = 1;
        }
        else
        {
            for (i = 0; i < 2; i++) //replace ph by pstar
            {
                ph.position[i] = pstar.position[i];
            }
        }
    }

    if (b==1) //Decision 'c'
    {
        if (ystar < y_h )
        {
            for (i = 0; i < 2; i++) //replace ph by pstar
            {
                ph.position[i] = pstar.position[i];
            }
        }
        //calculate pstarstar2
        for (i = 0; i < 2; i++)
        {
            pstarstar.position[i] = pstarstar2f(ph.position, pbar.position, i);
        }
        ystarstar = rosenbrockf(pstarstar.position);
        c = 1;
    }

    if (c == 1) //Decision 'e'
    {
        if (ystarstar > y_h)
        {
            for (i = 0; i < 3; i++) //replace all p_i's with (pi+pl)/2
            {
                for (i2 = 0; i2 < 2; i2++) // for pl // redundant
                {
                    pl.position[i2] = pbarf(pl.position, pl.position, i2);
                }
                for (i2 = 0; i2 < 2; i2++) //for pm
                {
                    pm.position[i2] = pbarf(pl.position, pm.position, i2);
                }
                for (i2 = 0; i2 < 2; i2++) //for ph
                {
                    ph.position[i2] = pbarf(pl.position, ph.position, i2);
                }
            }    
        }
        else
        {
            for (i = 0; i < 2; i++) //replace ph by pstarstar
            {
                ph.position[i] = pstarstar.position[i];
            }
        }
        
    }

    y_l = rosenbrockf(pl.position);
    y_m = rosenbrockf(pm.position);
    y_h = rosenbrockf(ph.position); 

    ybar = (y_l + y_m + y_h) / 3;

    for (i = 0; i < 2; i++) //redefines p1, p2, p3
    {
        p1.position[i] = pl.position[i];
    }
    for (i = 0; i < 2; i++)
    {
        p2.position[i] = pm.position[i];
    }
    for (i = 0; i < 2; i++)
    {
        p3.position[i] = ph.position[i];

    }
    
    } while(count <1000 && standevf(y_l, y_m, y_h, ybar)); //if true, will continue iterating
    
    printf("\ncount = %lf\n", *pcount);
    
    printf("\n\n");
    printf("pl position = %lf, ", pl.position[0]);
    printf("%lf\n", pl.position[1]);

    printf("\n\n");
    printf("pm position = %lf, ", pm.position[0]);
    printf("%lf\n", pm.position[1]);

    printf("\n\n");
    printf("ph position = %lf, ", ph.position[0]);
    printf("%lf\n\n", ph.position[1]);    

    printf("Elevation of minimum = %lf\n\n", rosenbrockf(ph.position));
    
    return 0;
}

int main()
{
    algorithm();  
    write_file("data.csv");

    return 0;
}