#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "LDPC_mod2.h"



//This file is for making parity-check matrix;
//These functions are all from folder "construct_LDPC", which is from Xuwei and NingdeXie.
//For original functions, please refer to "construct_LDPC".

int LDPC_mod2::demo_make_H_matrix()
{
	main_construct_H_matrix();

return 1;
	pchk_file="./H/LDPC_2x24ex382w4_20100830pm_2";
	gen_file="./H/LDPC_2x24ex382w4_20100830pm_2.g";	
	make_gen_main(pchk_file,gen_file,"sparse",0,0,0,0);

}

// To generate irregular LDPC according to given weight distributions.
int LDPC_mod2::main_construct_H_matrix()
{//this file is modified a little by dgq on 20091130.

	//H: (ex_factor*row_base_H) by (ex_factor*column_base_H);
	//H: column weight: sub_weight * row_base_H;

   // for quasi-cyclic LDPC

   char filename[200];
   int nf = 1;   //the number of files/ldpc codes.

	int index, i, j;
   WD_vector *row_dt, *col_dt;
 

   //the element of base matrix is a small circulant matrix.


   //input:
   int ex_factor = 384;//86; //the size of submatrix/basic square matrix.
   int sub_weight=2; //weight of circulant matrix;
   int row_base_H=2; //row weight of base matrix;
   int column_base_H=24; //column weight of base matrix;


   srand( (unsigned)time(NULL) ); // Ramdomize seed
   printf("generate %d quasi-cyclic LDPC codes: \n", nf);
   sprintf(filename, "LDPC_2x24ex382w4_20100830pm_2.1");
//   pchk_file="./H/LDPC_2x24ex382w4_20100830pm_2";

//	gen_quasi_ldpc_files(filename, nf, 3, 51, 1, ex_factor);
   gen_quasi_ldpc_files(filename, nf, row_base_H, column_base_H, sub_weight, ex_factor); //H: (ex_factor*4) by (ex_factor*68);

	
	printf("\b begin to free memory.\b");
	return 0;

	free(col_dt->wd);
	free(row_dt->wd);
	free(row_dt);
	free(col_dt);

	return 0;
}


void LDPC_mod2::gen_quasi_ldpc_files(char *codefile,  /* File name for the LDPC codes to be constructed */
								  int filenum,     /* Number of LDPC codes to be constructed */
                          int b_rows,      /* Number of rows in base matrix */ //=>number of (basic matrix)/submatrix in rows of H.
						  // base matrix is the mother matrix, with element as a small circulant matrix.
								  int b_cols,      /* Number of columns in base matrix */  //=>number of basic matrix in columns of H.
								  int sub_w,       /* Weight for the circulant */
								  int ex_factor		//=> the size of submatrix.
								  )
{
   char file[200];
   mod2sparse *exMatrix;
   
   int i, file_counter;
   FILE *f;   
	char postfix[5];
      
   file_counter = 0;
   printf("%d LDPC Codes to be generated ... \n", filenum);
   while (file_counter < filenum)
   {
      exMatrix = mod2sparse_allocate(b_rows * ex_factor, b_cols * ex_factor);
		if ( quasi_cyclic(exMatrix, b_rows, b_cols, sub_w, ex_factor) )
      {           
         sprintf( postfix, "%d", file_counter+1 ); 
         sprintf(file, "%s", codefile);
         sprintf( file + strlen(codefile), "%s", postfix);

         // Write expanded matrix to a file

         f = open_file_std(file,"wb");
         if (f==NULL) 
         { 
            fprintf(stderr,"Can't create parity check file for expanded matrix: %s\n", file);
             exit(1);
         }

         intio_write(f,('P'<<8)+0x80);

         if (ferror(f) || !mod2sparse_write(f,exMatrix) || fclose(f)!=0)
         { 
            fprintf(stderr,"Error writing to expanded matrix file %s\n",file);
            exit(1);
         } 
			
			file_counter++;
      }
		mod2sparse_free(exMatrix);
		
   }
}


int  LDPC_mod2::quasi_cyclic ( 
					mod2sparse *exPCH,      /* place to store the expanded ldpc matrix */
                    int rownum,             /* Number of rows in base matrix */
                    int colnum,             /* Number of columns in base matrix */
					int sub_w,              /* weight of each circulant submatrix */
                    int ex_factor           /* expansion factor */
						 ) 

{//the core function to generate QC-ldpc code:
   int i, j, k, k1, p1;
	int cnt;
   
   tanner_graph *ex_g;// *ex_g_aux; // tanner graph for the expanded matrix
   int shift;	
   int cycle_constraint;

   //int shift_value[CNU_NUM][VNU_NUM][BLK_WEIGHT];
   int ***shift_value; //rownum*colnum*sub_w;
   shift_value=(int ***)calloc(rownum,sizeof(int **));
   for(i=0;i<rownum;i++)
   {
	   shift_value[i] = (int **) calloc(colnum,sizeof(int *));
	   for(j=0;j<colnum;j++)
		   shift_value[i][j] = (int *)calloc(sub_w,sizeof(int));
   }

   // Initialize tanner graph for Breath-first search in checking cycle and cycle degree
   ex_g = (tanner_graph *)calloc(1, sizeof(tanner_graph));
   ex_g->bit_node = (node_entry *)calloc(colnum*ex_factor, sizeof(node_entry));
   ex_g->check_node = (node_entry *)calloc(rownum*ex_factor, sizeof(node_entry));
   ex_g->matrix = exPCH;

   // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

   srand( (unsigned)time(NULL) ); // Ramdomize seed
	

   for(i = 0; i < rownum; i++)
   {
		printf("\n\n i = %d; \n", i);
      for(j = 0; j < colnum; j++)
      {         
			//printf("\n j = %d; \n", j);
         for(k = 0; k < sub_w; k++)
         {  
				cnt = 0;            
retry:      shift = rand() % ex_factor;				
				
				for(k1 = 0; k1 < k; k1++)
				{
					if (shift == shift_value[i][j][k1])
						goto retry;
				}
				shift_value[i][j][k] = shift;
            insert_sub_matrix( exPCH, i * ex_factor,  j * ex_factor, ex_factor, shift );
				                  
				cycle_constraint = 0;

				if(i==3 && j==37)
				{
					i=i;
				}
				printf("i=%d;j=%d\n",i,j);
				
				for(p1 = 0; p1 < ex_factor; p1++)
            {
					if ( check_cycle(ex_g, i * ex_factor + p1, 4) )                     
               {
                   cycle_constraint = 1; //short cycle is found.
                    break;
                }                     
			}

            cnt++;
			if ( cycle_constraint )
            {
					remove_sub_matrix( exPCH, i * ex_factor, j * ex_factor, ex_factor, shift );
					if (cnt > 20*ex_factor) //if have tried too many times.
						return 0;
					else
						goto retry; //
				}			
				//printf(" k = %d; %d ", k, shift);	
			}
		}            
   }

   // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	

	for(i = 0; i < rownum; i++)
	{
      for(j = 0; j < colnum; j++)
			for(k = 0; k < sub_w; k++)
				init_addr[i*colnum*sub_w + j*sub_w + k] = shift_value[i][j][k];
	}

   free(ex_g->bit_node);
   free(ex_g->check_node);
   free(ex_g);
   
 
  for(i=0;i<rownum;i++)
  {
	  for(j=0;j<colnum;j++)
		free(shift_value[i][j]);

	  free(shift_value[i]);
  }
  free(shift_value);
 
  
   return 1;

}


void LDPC_mod2::insert_sub_matrix( mod2sparse *mainMatrix,
                        int startRow,
                        int startCol,
                        int subDim,
                        int shift )

{
   int i, rowPos, colPos;
   for( i = 0; i < subDim; i++)
   {
      rowPos = startRow + i;
      colPos = startCol + (i + shift)%subDim;
      mod2sparse_insert(mainMatrix, rowPos, colPos);
   }
      
}


//-----------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------
void LDPC_mod2::remove_sub_matrix( mod2sparse *mainMatrix,
                        int startRow,
                        int startCol,
                        int subDim,
                        int shift )
{
   int i, rowPos, colPos;
   for( i = 0; i < subDim; i++)
   {
      rowPos = startRow + i;
      colPos = startCol + (i + shift)%subDim;
      mod2sparse_delete(mainMatrix, mod2sparse_find(mainMatrix, rowPos, colPos));
   }
}

//----------
//----------------------------------------------------------------------------
// assume the root node is always check node 
int LDPC_mod2::check_cycle(tanner_graph *g, int root_index, int length)
{
   int M, N, i;
   int current_index, tmp_index;
   node_entry *current_node, *tmp_node;
   node_entry p_node1, p_node2;  // used for tracing back
   mod2entry *e;
   Queue Q;

   M = mod2sparse_rows(g->matrix);
   N = mod2sparse_cols(g->matrix);
     
   Q = CreateQueue( M + N );
   
   // Initialize check nodes 
   for (i = 0; i < M; i++)
   {
      g->check_node[i].index = i;
      g->check_node[i].color = white;
      g->check_node[i].parent = NULL;
      g->check_node[i].depth = -1;
   }

   // Initialize bit nodes
   for (i = 0; i < N; i++)
   {
      g->bit_node[i].index = i + M;
      g->bit_node[i].color = white;
      g->bit_node[i].parent = NULL;
      g->bit_node[i].depth = -1; 
   }

   g->check_node[root_index].color = gray;
   g->check_node[root_index].depth = 0;
   Enqueue(root_index, Q);

   while( !IsEmpty( Q ) )
   {
      current_index = FrontAndDequeue( Q );

      if (current_index < M) // Current node is check node
      {
         current_node = &g->check_node[current_index];
         // get the first bit node connected to this check node
         e = mod2sparse_first_in_row(g->matrix, current_index);          
         
         while (!mod2sparse_at_end(e)) 
	      {
            //retrieve the bit node index
            tmp_index = mod2sparse_col(e);
            tmp_node = &g->bit_node[tmp_index];

            // If the bit node is not the parent of the check node, then don't try to
            // to visit it. If not, there are two cases to analyse: if it is in white 
            // color then visit it; if it's gray and the ancestors right under the root
            // are different, then the root node is in a cycle. 
            if(current_node->parent == NULL || tmp_index != current_node->parent->index)
            {
               switch(tmp_node->color)
               {
                  case white :               
                     tmp_node->color = gray;
                     tmp_node->parent = current_node;
                     tmp_node->depth = 1 + current_node->depth;
                     Enqueue( tmp_node->index, Q);
                     break;
                              
                  case gray :               
                     if (2*tmp_node->depth <= length)
                     {
                        p_node1 = trace_back(*current_node, current_node->depth - 1);
                        p_node2 = trace_back(*tmp_node, tmp_node->depth - 1);
                        if (p_node1.index != p_node2.index)
                        {
                           DisposeQueue(Q);
                           return 1; // find cycle with given length passing root
                        }
                     } 
                     // esle the cycle does not pass the root
                     break;

                  case black:
                     break;  // ignore, the edge is already used

                  default: 
                     fprintf(stderr, "Invalid node color!");
               }               
            }

            // get next bit node connected to this check node
            e = mod2sparse_next_in_row(e);

         } 
         
         current_node->color = black;
      }

      if (current_index >= M) // Current node is bit node 
      {
         current_index -=  M;
         current_node = &g->bit_node[current_index];

         // get the first check node connected to this bit node
         e = mod2sparse_first_in_col(g->matrix, current_index);
         
         while (!mod2sparse_at_end(e))
	      {	        
            // retrieve the check node index
            tmp_index = mod2sparse_row(e);
            tmp_node = &g->check_node[tmp_index];

            // if this check node is not the parent of the bit node, then there are
            // two cases to analyse: 
            // if it is in white color then visit it; 
            // if it's already in gray color, then we assert that it has formed a cycle.  
            if( current_node->parent == NULL || tmp_index != current_node->parent->index)
            {
               switch (tmp_node->color)
               {
               case white :
                  tmp_node->color = gray;
                  tmp_node->parent = current_node;
                  tmp_node->depth = 1 + current_node->depth;
                  Enqueue(tmp_node->index, Q);
                  break;

               case gray :               
                  if (2*tmp_node->depth <= length)
                  {
                     p_node1 = trace_back(*current_node, current_node->depth - 1);
                     p_node2 = trace_back(*tmp_node, tmp_node->depth - 1);
                     if (p_node1.index != p_node2.index)
                     {
                        DisposeQueue(Q);
                        return 1; // find cycle with given length passing root
                     }
                     // esle the cycle does not pass the root
                  }
                  break; 
               } 

            } 

            // get next check node connected to this bit node
            e = mod2sparse_next_in_col(e);
         }           
         current_node->color = black;
      
      }

      if ( 2 * current_node->depth > length )
      {
         DisposeQueue(Q);
         return 0;
      }
   
   }

   DisposeQueue(Q);
   return 0;
}


//-----------------------------------------------------------------


int LDPC_mod2::IsEmpty( Queue Q )
{
   return Q->Size == 0;
}

int LDPC_mod2::IsFull( Queue Q )
{
   return Q->Size == Q->Capacity;
}


Queue LDPC_mod2::CreateQueue( int MaxElements )
{
   Queue Q;

   if( MaxElements < MinQueueSize )
   {
	   printf( "Queue size is too small" );
	   exit(0);
   }

   //?
   Q = (QueueRecord *) calloc(1, sizeof( struct QueueRecord ) );
   if( Q == NULL )
   {
	   printf( "Out of space!!!" );
	   exit(0);
   }

   //?
   Q->Array = (ElementType*)calloc(MaxElements ,sizeof( ElementType ) );
   if( Q->Array == NULL )
   {
	   printf( "Out of space!!!" );
		exit(0);
   }
   Q->Capacity = MaxElements;
   MakeEmpty( Q );

   return Q;
}

void LDPC_mod2::MakeEmpty( Queue Q )
{
   Q->Size = 0;
   Q->Front = 1;
   Q->Rear = 0;
}


void LDPC_mod2::DisposeQueue( Queue Q )
{
   if( Q != NULL )
   {
       free( Q->Array );
       free( Q );
   }
}

 int LDPC_mod2::Succ( int Value, Queue Q )
{
   if( ++Value == Q->Capacity )
       Value = 0;
   return Value;
}

void LDPC_mod2::Enqueue( ElementType X, Queue Q )
{
   if( IsFull( Q ) )
   {
	   printf( "Full queue" );
	   exit(0);
   }
   else
   {
       Q->Size++;
       Q->Rear = Succ( Q->Rear, Q );
       Q->Array[ Q->Rear ] = X;
   }
}

ElementType LDPC_mod2::Front( Queue Q )
{
   if( !IsEmpty( Q ) )
       return Q->Array[ Q->Front ];

   printf( "Empty queue" );
   exit(0);
   
   return 0;  /* Return value used to avoid warning */
}

void LDPC_mod2::Dequeue( Queue Q )
{
   if( IsEmpty( Q ) )
   {
       printf( "Empty queue" );
	   exit(0);
   }
   else
   {
       Q->Size--;
       Q->Front = Succ( Q->Front, Q );
   }
}

ElementType LDPC_mod2::FrontAndDequeue( Queue Q )
{
   ElementType X = 0;

   if( IsEmpty( Q ) )
   {
       Error( "Empty queue" );
   }
   else
   {
       Q->Size--;
       X = Q->Array[ Q->Front ];
       Q->Front = Succ( Q->Front, Q );
   }
   return X;
}



/* OPEN A FILE THAT MIGHT BE STANDARD INPUT OR OUTPUT.  If the file name
   given is "-", this procedure just returns stdin or stdout, depending on
   whether the mode is for reading or writing.  Otherwise, fopen is called. 
*/

FILE * LDPC_mod2::open_file_std ( char *fname,	/* Name of file to open, or "-" for stdin/stdout */
                      char *mode	 )/* Mode for opening: eg, "r" or "w" */
{ 
   if (strcmp(fname,"-") == 0)
   { 
      switch (mode[0])
      { 
         case 'r': 
            return stdin;
         case 'w': 
            return stdout;
         default:  
            { 
               fprintf(stderr,"Bad mode passed to open_file_std: %s\n",mode);
               exit(1);
            }
      }
   }
   else
      return fopen(fname,mode);
}


node_entry LDPC_mod2::trace_back(node_entry u, int n)
{
   int i;
   for(i = 0; i < n; i++)
   {
      u = *u.parent;
   }
   return u;
}
