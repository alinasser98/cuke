import sys
import os

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from core.ir import *
from codegen.cpu import *

def PrintCCode(ir):
	code = ''
	for d in ir:
		if d:
			code += to_string(d)
	print(code)


def Loop0():
    ir = []

    N = Scalar('int', 'N')
    M = Scalar('int', 'M')
    L = Scalar('int', 'L')
    A = Ndarray('int', (N, M, L), 'A')
    B = Ndarray('int', (N, M, L), 'B')

    loopi = Loop(0, N, 1, [])
    loopj = Loop(0, M, 1, [])
    loopk = Loop(0, L, 1, [])

    loopi.body.append(loopj)
    loopj.body.append(loopk)

    lhs1 = Index(Index(Index(A, Expr(loopi.iterate, 1, '+')), loopj.iterate), loopk.iterate)
    lhs2 = Index(Index(Index(B, Expr(loopi.iterate, 1, '+')), Expr(loopj.iterate, 2, '+')), Expr(loopk.iterate, 1, '-'))
    rhs1 = Index(Index(Index(B, Expr(loopi.iterate, 1, '+')), loopj.iterate), Expr(loopk.iterate, 1, '-'))
    rhs2 = Index(Index(Index(A, loopi.iterate), loopj.iterate),  Expr(loopk.iterate, 1, '+'))
    rhs3 = Index(Index(Index(B, loopi.iterate), Expr(loopj.iterate, 2, '+')),  loopk.iterate)
	
    # body = Assignment(lhs, Expr(rhs1, rhs2, '+'))
    loopk.body.extend([Assignment(lhs1, Expr(rhs1, 2, '+')), Assignment(lhs2, Expr(rhs2, rhs3, '+'))])

    ir.extend([Decl(L)])
    ir.extend([Decl(M)])
    ir.extend([Decl(N)])
    ir.extend([Decl(A)])
    ir.extend([loopi])

    return ir


# for ( k = 0; k < L ; ++k ){
# 	for ( j = 0; j < M; ++ j ){
# 		for ( i = 0; i < N; ++ i ){
# 			a[i+1] [j+1] [k] = a [i] [j] [k] + a [i] [j + 1] [k + 1] ;
# 		}
# 	}
# }

#Distance Vector: 
#[1, 1, 0] :  a[i+1] [j+1] [k] and a [i] [j] [k]
#[1, 0, -1] : a[i+1] [j+1] [k] and  a [i] [j + 1] [k + 1] 

#Direction Vector:
#[<, <, =]
#[<, =, >]

def Loop1():
    ir = []

    L = Scalar('int', 'L')
    M = Scalar('int', 'M')
    N = Scalar('int', 'N')
    A = Ndarray('int', (N, M, L), 'A')

    loopk = Loop(0, L, 1, [])
    loopj = Loop(0, M, 1, [])
    loopi = Loop(0, N, 1, [])
    loopk.body.append(loopj)
    loopj.body.append(loopi)

    lhs = Index(Index(Index(A, Expr(loopi.iterate, 1, '+')), Expr(loopj.iterate, 1, '+')), loopk.iterate)
    rhs1 = Index(Index(Index(A, loopi.iterate), loopj.iterate), loopk.iterate)
    rhs2 = Index(Index(Index(A, loopi.iterate), Expr(loopj.iterate, 1, '+')),  Expr(loopk.iterate, 1, '+'))
	
    body = Assignment(lhs, Expr(rhs1, rhs2, '+'))
    loopi.body.append(body)

    ir.extend([Decl(L)])
    ir.extend([Decl(M)])
    ir.extend([Decl(N)])
    ir.extend([Decl(A)])
    ir.extend([loopk])

    return ir


# for ( i = 0; i < N ; ++i ){
# 	for ( j = 0; j < N; ++ j ){
# 			a[i][j] = a[i+1][j-1];
# 	}
# }

#Distance Vector: 
#[-1, 1]

#Direction Vector:
#[<, >]

def Loop2():
    ir = []

    N = Scalar('int', 'N')
    A = Ndarray('int', (N, N), 'A')

    loopi = Loop(0, N, 1, [])
    loopj = Loop(0, N, 1, [])

    loopi.body.append(loopj)

    lhs = Index(Index(A, loopi.iterate), loopj.iterate)
    rhs = Index(Index(A, Expr(loopi.iterate, 1, '+')), Expr(loopj.iterate, 1, '-'))
	
    loopj.body.append(Assignment(lhs, rhs))

    ir.extend([Decl(N)])
    ir.extend([Decl(A)])
    ir.extend([loopi])

    return ir

######################################################## Steps
# 1. Identify the loop body
# the output of first step
# A[_l0 + 1][_l1][_l2] = B[_l0 + 1][_l1][_l2 - 1] + 2;
# B[_l0 + 1][_l1 + 2][_l2 - 1] = A[_l0][_l1][_l2 + 1] + B[_l0][_l1 + 2][_l2];


# 2. Get the index statement of the loop body
# write ststement array: A[_l0 + 1][_l1][_l2]B[_l0 + 1][_l1 + 2][_l2 - 1]
# read statement array: B[_l0 + 1][_l1][_l2 - 1]A[_l0][_l1][_l2 + 1]B[_l0][_l1 + 2][_l2]

# 3. Group the index statement by there names.
# Two dict:
#     Write dicts: {A : [A[_l10 + 1][_l1][_l2]], 'B' : [B[_l0 + 1][_l1 + 2][_l2 - 1]]}
#     Read dicts:  {A : [A[_l0][_l1][_l2 + 1]], 'B' : [B[_l0 + 1][_l1][_l2 - 1], B[_l0][_l1 + 2][_l2]]}

# 4. Compute the direction vector
# [Write expression, Read expression]
# [A[_l10 + 1][_l1][_l2], A[_l0][_l1][_l2 + 1]]
# [B[_l0 + 1][_l1 + 2][_l2 - 1], B[_l0 + 1][_l1][_l2 - 1]]
# [B[_l0 + 1][_l1 + 2][_l2 - 1], B[_l0][_l1 + 2][_l2]]

# Distance vector:
# [1, 0, -1]
# [0, 2,  0]
# [1, 0, -1]

# Direction vector
# [<, =, >]
# [=, <, =]
# [<, =, >]

# 5. Safety checking based on the direction vector we get:
# Exchange [0, 1]
# [=, <, >]
# [<, =, =]
# [=, <, >]

# Exchange [0, 2]
# [>, =, <] The first 

######################################################## GetIndex()
def GetIndex(statement, is_write, write_expr = [], read_expr = []):
    if type(statement) == Ndarray or type(statement) == Index:
        if is_write:
            write_expr.append(statement)
        else:
            read_expr.append(statement)
            
    if type(statement) == Assignment:
        GetIndex(statement.lhs, True, write_expr, read_expr)
        GetIndex(statement.rhs, False, write_expr, read_expr)
    elif type(statement) == Expr:
        GetIndex(statement.left, is_write, write_expr, read_expr)
        GetIndex(statement.right, is_write, write_expr, read_expr)
    else:
        return
######################################################## FindBody() Part 1
def FindBody(nested_loop):
    if not type(nested_loop) == Loop:
        return nested_loop
    if type(nested_loop.body[0]) == Loop:
        return FindBody(nested_loop.body[0])
    else:
        return nested_loop.body
######################################################## InterchangeLoop()
def InterchangeLoop(ir, loop_idx=[]):
    
    ir_res = []
    write_expr = []
    read_expr = []
    #######################################
    for ir_item in ir:
        if type(ir_item)==Loop:
            body = FindBody(ir_item)
            for body_item in body:
                GetIndex(body_item, False, write_expr, read_expr)
    
    
    #######################################     
    print("<===== Loop body we got! =====>")    
    PrintCCode(body)
    
    print("<===== Read/Write Statements =====>")
    # Index data structure
    PrintCCode(write_expr)  
    PrintCCode(read_expr)


if __name__ == "__main__":
    loop0_ir = Loop0()
    loop1_ir = Loop1()
    loop2_ir = Loop2()
    PrintCCode(loop0_ir)

    optimized_loop1_ir = InterchangeLoop(loop0_ir, [0, 1])
    # optimized_loop1_ir = InterchangeLoop(loop1_ir, [1, 2])
    # optimized_loop2_ir = InterchangeLoop(loop2_ir, [0, 1])

    # optimized_ir = LoopInterchange(ir)
    # print("Loop after interchange:")
    # PrintCCode(optimized_ir)
