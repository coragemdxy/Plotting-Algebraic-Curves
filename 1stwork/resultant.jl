using Nemo

# ============================================================
# 构造 F 和 G 关于主变量的 Sylvester 矩阵
#
# 输入：
#   F, G：同一个一元多项式环中的两个多项式
#
# 输出：
#   Sylvester 矩阵
# ============================================================

function sylvester_matrix(F, G)
    # 确保 F 和 G 属于同一个多项式环
    parent(F) == parent(G) ||
        throw(ArgumentError("F 和 G 必须属于同一个一元多项式环"))

    # 零多项式没有通常意义下的次数
    iszero(F) && throw(ArgumentError("F 不能是零多项式"))
    iszero(G) && throw(ArgumentError("G 不能是零多项式"))

    # F 和 G 的次数
    m = degree(F)
    n = degree(G)

    # 系数所在的环
    R = base_ring(parent(F))

    # Sylvester 矩阵的阶数
    N = m + n

    # 创建 N × N 的零矩阵
    S = zero_matrix(R, N, N)

    # --------------------------------------------------------
    # 前 n 行放置 F 的系数
    #
    # F = a_m t^m + a_{m-1} t^{m-1} + ... + a_0
    #
    # 每一行向右移动一列
    # --------------------------------------------------------

    for row in 1:n
        for j in 0:m
            S[row, row + j] = coeff(F, m - j)
        end
    end

    # --------------------------------------------------------
    # 后 m 行放置 G 的系数
    #
    # G = b_n t^n + b_{n-1} t^{n-1} + ... + b_0
    #
    # 每一行向右移动一列
    # --------------------------------------------------------

    for row in 1:m
        for j in 0:n
            S[n + row, row + j] = coeff(G, n - j)
        end
    end

    return S
end


# ============================================================
# 计算 Sylvester 结式
#
# 不调用 resultant。
# 结式按照定义等于 Sylvester 矩阵的行列式。
# ============================================================

function sylvester_resultant(F, G)
    S = sylvester_matrix(F, G)
    return det(S)
end


# ============================================================
# 测试题目中的例子
#
# res((t^2 + 1)x - t, (t^2 + 1)y - 2, t)
# ============================================================

function test_example()
    # 创建系数环 Q[x,y]
    Rxy, variables = polynomial_ring(QQ, ["x", "y"])
    x, y = variables

    # 在 Q[x,y] 上创建关于 t 的一元多项式环
    Rt, t = polynomial_ring(Rxy, "t")

    # 定义两个关于 t 的一元多项式
    F = (t^2 + 1) * x - t
    G = (t^2 + 1) * y - 2

    println("F = ")
    println(F)
    println()

    println("G = ")
    println(G)
    println()

    S = sylvester_matrix(F, G)

    println("Sylvester 矩阵为：")
    println(S)
    println()

    result = sylvester_resultant(F, G)

    println("计算得到的结式为：")
    println(result)
    println()

    expected = 4*x^2 + y^2 - 2*y

    println("题目给出的结果为：")
    println(expected)
    println()

    println("程序验证是否正确：")
    println(result == expected)

    # --------------------------------------------------------
    # 第一步：建立系数环 Q[x,y,z,t]
    #
    # 为了先关于 s 做结式，
    # 把 x,y,z,t 看作系数，s 看作主变量
    # --------------------------------------------------------

    Rxyzt, vars = polynomial_ring(QQ, ["x", "y", "z", "t"])
    x, y, z, t = vars

    Rs, s = polynomial_ring(Rxyzt, "s")

    F2 = (s^2+1)y - (s^2-1)*t^2
    F3 = (s^2+1)z -2s*t^2

    println("F2 = ")
    println(F2)
    println()

    println("F3 = ")
    println(F3)
    println()

    # --------------------------------------------------------
    # 第二步：消去 s
    # --------------------------------------------------------

    R_s = sylvester_resultant(F2, F3)

    println("关于 s 的结式为：")
    println(R_s)
    println()

    # 理论上：
    #
    # R_s = 4*t^4*(y^2 + z^2 - t^4)
    #
    # 去掉无关因子 4*t^4，得到真正需要的关系
    H_base = y^2 + z^2 - t^4

    println("去掉无关因子后得到：")
    println(H_base)
    println()

    # --------------------------------------------------------
    # 第三步：重新建立关于 t 的一元多项式环
    #
    # 此时 x,y,z 是系数，t 是主变量
    # --------------------------------------------------------

    Rxyz, vars_xyz = polynomial_ring(QQ, ["x", "y", "z"])
    x2, y2, z2 = vars_xyz

    Rt, t2 = polynomial_ring(Rxyz, "t")

    F1 = 2*x2 - t2^3
    H  = y2^2 + z2^2 - t2^4

    println("F1 = ")
    println(F1)
    println()

    println("H = ")
    println(H)
    println()

    # --------------------------------------------------------
    # 第四步：消去 t
    # --------------------------------------------------------

    implicit_polynomial = sylvester_resultant(F1, H)

    println("最终隐式多项式为：")
    println(implicit_polynomial)
    println()

    expected = (y2^2 + z2^2)^3 - 16*x2^4

    println("整理后的标准形式为：")
    println(expected)
    println()

    println("程序结果是否正确：")
    println(implicit_polynomial == expected)
end


# 只有直接运行本文件时才执行测试
if abspath(PROGRAM_FILE) == @__FILE__
    test_example()
end