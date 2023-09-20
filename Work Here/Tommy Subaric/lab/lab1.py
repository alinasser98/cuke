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
 
def PrintCGroupCode(index_groups):
    for name, expressions in index_groups.items():
        for expr in expressions:
            print(to_string(expr))

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
#   for ( j = 0; j < M; ++ j ){
#       for ( i = 0; i < N; ++ i ){
#           a[i+1] [j+1] [k] = a [i] [j] [k] + a [i] [j + 1] [k + 1] ;
#       }
#   }
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
#   for ( j = 0; j < N; ++ j ){
#           a[i][j] = a[i+1][j-1];
#   }
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

# 1. Identify loop body
# A[_l0 + 1][_l1][_l2] = B[_l0 + 1][_l1][_l2 - 1] + 2;
# B[_l0 + 1][_l1 + 2][_l2 - 1] = A[_l0][_l1][_l2 + 1] + B[_l0][_l1 + 2][_l2];

def FindBody(nested_loop):
    if type(nested_loop) != Loop:
        return nested_loop.body
    if type(nested_loop.body[0]) == Loop:
        return FindBody(nested_loop.body[0])
    else:
        return nested_loop.body

# 2. get the index statement of the loop body
# Write statement array: [A[_l0 + 1][_l1][_l2], B[_l0 + 1][_l1 + 2][_l2 - 1]]
# Read statement array: [B[_l0 + 1][_l1][_l2 - 1], A[_l0][_l1][_l2 + 1], B[_l0][_l1 + 2][_l2]]

def GetIndex(statement, is_write, write_expr = None, read_expr = None):
    if read_expr is None:
        read_expr = list()
    if write_expr is None:
        write_expr = list()
    if type(statement)==Ndarray or type(statement)== Index:
        if is_write:
            write_expr.append(statement)
        else:
            read_expr.append(statement)
    if type(statement)==Assignment:
        GetIndex(statement.lhs, True, write_expr, read_expr)
        GetIndex(statement.rhs, False, write_expr, read_expr)
    elif type(statement)==Expr:
        GetIndex(statement.left, is_write, write_expr, read_expr)
        GetIndex(statement.right, is_write, write_expr, read_expr)
    else:
        return

# 3. Group the index statement by their names. 
# Two dict:
# Write dict: {'A':[A[_l0 + 1][_l1][_l2], 'B': B[_l0 + 1][_l1][_l2 - 1]}
# Read dict: {'A': [A[_l0][_l1][_l2 + 1]], 'B': [B[_l0 + 1][_l1][_l2 - 1] , B[_l0][_l1 + 2][_l2]]}
    
def GroupIndicesByNames(write_expr, read_expr):
    def recursive_group(statement, write_dict, read_dict):
        if isinstance(statement, (Ndarray, Index)):
            name = statement.name  # Assuming there is a 'name' attribute in Index or Ndarray
            if 'A' not in write_dict:
                write_dict[name] = []
            if 'B' not in read_dict:
                read_dict[name] = []

            if statement in write_expr:
                write_dict[name].append(statement)
            elif statement in read_expr:
                read_dict[name].append(statement)
        elif isinstance(statement, Assignment):
            recursive_group(statement.lhs, write_dict, read_dict)
            recursive_group(statement.rhs, write_dict, read_dict)
        elif isinstance(statement, Expr):
            recursive_group(statement.left, write_dict, read_dict)
            recursive_group(statement.right, write_dict, read_dict)

    write_index_groups = {}
    read_index_groups = {}

    for statement in write_expr + read_expr:
        recursive_group(statement, write_index_groups, read_index_groups)

    return write_index_groups, read_index_groups

#4. Compute the direction vector
# [A[_l0 + 1][_l1][_l2], A[_l0][_l1][_l2 + 1]]
# [B[_l0 + 1][_l1][_l2 - 1], B[_l0 + 1][_l1][_l2 - 1]]
# [B[_l0 + 1][_l1][_l2 - 1], B[_l0][_l1 + 2][_l2]]

# Direction Vector
# [<, =, >]
# [=, <, =]
# [<, =, >]

def ComputeDirectionVector(distance_vector):
    direction_vector = []
    for distance in distance_vector:
        if distance < 0:
            direction_vector.append('<')
        elif distance == 0:
            direction_vector.append('=')
        else:
            direction_vector.append('>')
    return direction_vector

# Distance Vector:
# [1,0,-1]
# [0,2,0]
# [1,0,-1]

def ComputeDistanceVector(expr1, expr2):
    distance_vector = []
    # Assuming expr1 and expr2 have the same shape
    for dim1, dim2 in zip(expr1.indices, expr2.indices):
        distance = dim1 - dim2
        distance_vector.append(distance)
    return distance_vector

def InterchangeLoop(loop, loop_idx=[]):
    ir_res = []
    write_expr = []
    read_expr = []

    for ir_item in loop:
        if isinstance(ir_item, Loop):
            body = FindBody(ir_item)
            for body_item in body:
                GetIndex(body_item, False, write_expr, read_expr)
                
    # Print the expressions          
    print('Write Expression:')
    PrintCCode(write_expr)
    print('Read Expression:')
    PrintCCode(read_expr)          

    # Calculate the index groups once after collecting all relevant expressions
    write_index_groups, read_index_groups = GroupIndicesByNames(write_expr, read_expr)

    # Print the index groups
    print('Write Index Groups:')
    PrintCGroupCode(write_index_groups)
    print('Read Index Groups:')
    PrintCGroupCode(read_index_groups)
    
if __name__ == "__main__":
    loop0_ir = Loop0()
    loop1_ir = Loop1()
    loop2_ir = Loop2()
    PrintCCode(loop0_ir)

    optimized_loop1_ir = InterchangeLoop(loop0_ir, [0, 1])
    
    # PrintCCode(optimized_loop1_ir)
    # optimized_loop1_ir = InterchangeLoop(loop1_ir, [1, 2])
    # optimized_loop2_ir = InterchangeLoop(loop2_ir, [0, 1])

    # optimized_ir = LoopInterchange(ir)
    # print("Loop after interchange:")
    # PrintCCode(optimized_ir)

# Safety Checking based on the direction vector we get. 
# Exchange [0, 1]
# [=, <, >]
# [<, =, =]
# [=, <, >]
# it is safe

# Exchange [0, 2]
# [>, =, <]
# [=, <, =]
# [>, =, <]
# it is NOT safe