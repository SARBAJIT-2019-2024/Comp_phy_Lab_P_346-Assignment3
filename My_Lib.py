#         #   #     #         #            #     # # # #        # # # #      #     #
# #     # #    #   #          #            #     #     #        #     #       #   #
#  #   #  #     # #           #            #     #     #        #     #        # #
#   # #   #      #            #            #     # # # #        # # # #         #
#    #    #      #            #            #     #     #        #               #
#         #      #            #            #     #     #        #               #  
#         #      #   #######  # # # # #    #     # # # #   #    #               # 




## Custom library for importing functions
##################################################################################################
#Partial Pivot
def partial_pivoting(a,b,n):
    n=len(a)
    no_of_swaps=0
    for i in range(n-1):
        if abs(a[i][i]) == 0:
            for j in range(i+1,n):
                if abs(a[j][i]) > abs(a[i][i]):
                    a[j], a[i] = a[i], a[j]  # interchange ith and jth rows of matrix 'A'
                    b[j], b[i] = b[i], b[j]  # interchange ith and jth elements of vector 'b'
            no_of_swaps+=1
    return no_of_swaps        

###################################################################################################
#defining Gauss Jordan
def Gauss_jordan(list_C):
    a,b=Matrixmaker(list_C)
    n=len(a)
    
    partial_pivoting(a,b,n) #calling the partial pivot
    
    for i in range(n):
        pivot_element = a[i][i]
        
        
        if pivot_element==0:
            return (None,None)
        if type(b[i]) is list: 
            for l in range(n):
                if b[i][l] != 0:
                    b[i][l] = b[i][l]/pivot_element
            for r in range(i, n):
                a[i][r] = a[i][r]/pivot_element
            for k in range(n):
                if k != i and a[k][i] != 0:
                    balance_factor = a[k][i]
                    for j in range(i,n):
                        a[k][j] = a[k][j] - balance_factor*a[i][j]  
                    for l in range(n):
                        if b[i][l] != 0:
                            b[k][l] = b[k][l] - balance_factor*b[i][l]

                        
        else: 
            b[i] = b[i]/pivot_element
            for k in range(i, n):
                 a[i][k] = a[i][k]/pivot_element

            for k in range(n):
                if k != i and a[k][i] != 0:
                    balance_factor = a[k][i]
                    b[k] = b[k] - balance_factor*b[i]
                    for j in range(i, n):
                        a[k][j] = a[k][j] - balance_factor*a[i][j] 
    return(a,b)
###############################################################################################
# Defining matrix maker which makes a matrix C in to A and b
    
def Matrixmaker(list_C):
    list_A=[[0 for x in range(len(list_C))] for y in range(len(list_C))]
    for i in range(len(list_C)):
            for j in range(len(list_C)):
                list_A[i][j]=list_C[i][j] 
    if len(list_C[0])==len(list_C)+1: 
        list_B=[0 for x in range(len(list_C))] 
        for i in range(len(list_C)):
            list_B.append(0)
        
        for i in range(len(list_C)):
            list_B[i]=list_C[i][len(list_C)]
    else:
        list_B=[[0 for x in range(len(list_C))] for y in range(len(list_C))] 
        for i in range(len(list_C)):

            list_B[i][i]=1       
    return(list_A,list_B)

#####################################################################################################
 #Defining matrix multiplication

def matrix_mul(a,b):
    AB_=[[0 for x in range(len(a))] for y in range(len(a))] 
        
    for i in range(len(a)):
        for j in range(len(b[2])):#expanding along the coloumn 2 of B matrix
            for k in range(len(b)):
                AB_[i][j] += a[i][k]*b[k][j] 
    return AB_
#####################################################################################################    

# Defining function for determinant calculation

def determinant_calc(a):           
    
        no_of_swaps = 0
        for i in range(len(a)): 
            if abs(a[i][i])== 0:
                for j in range(i+1 , len(a)): 
                    if  abs(a[i][i]) < abs(a[j][i]): 
                        for k in range(i , len(a)):
                            a[i][k], a[j][k] = a[j][k], a[i][k]         #swapping the rows of the matrix.
                no_of_swaps += 1  
        for i in range(len(a)):
            pivot_element = a[i][i]
            if pivot_element==0:
                return (0)
            for k in range(len(a)):
                if k != i and a[k][i] != 0:
                    balance_factor = a[k][i]
                    for j in range(i,len(a)):
                        a[k][j] = a[k][j] - balance_factor*(a[i][j]/pivot_element)
        det=1
        for i in range(len(a)):
            det*=a[i][i]
        det*=(-1)**(no_of_swaps)    
                
            
        return(det)     
                
###########################################################################################################
##############  
    












