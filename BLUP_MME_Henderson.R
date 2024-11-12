####BLUP_MME_Henderson####

#Função para calcular BLUP e BLUE - MLE
bb_MLE <- function(y, X, Z, sigma2_e, sigma2_u) {
  #Dimensões das matrizes
  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(Z)
  
  # Matriz de variância do resíduo
  R <- diag(sigma2_e, n, n)
  
  # Matriz de variância de efeitos aleatórios
  G <- diag(sigma2_u, q, q)
  
  # Covariance matrix for the vector of observation y
  V = (Z %*% G %*% t(Z)) + R
  
  # BLUE
  b = (solve(t(X) %*% solve(V) %*% X)) %*% (t(X) %*% solve(V) %*% y)
  
  # BLUP
  u = G %*% t(Z) %*% solve(V) %*% (y - (X %*% b))
  
  # Retornando os valores
  list(BLUE = b, BLUP = u)
}

#Função para calcular BLUP e BLUE - Henderson
bb_Hend1 <- function(y, X, Z, sigma2_e, sigma2_u) {
  #Dimensões das matrizes
  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(Z)
  
  # Matriz de variância do resíduo
  R <- diag(sigma2_e, n, n)
  
  # Matriz de variância de efeitos aleatórios
  G <- diag(sigma2_u, q, q)
  
  # Matrix parameters
  h_a = t(X) %*% solve(R) %*% X
  
  h_b = t(X) %*% solve(R) %*% Z
  
  h_c = t(Z) %*% solve(R) %*% X
  
  h_d = (t(Z) %*% solve(R) %*% Z) + solve(G)
  
  up = cbind(h_a, h_b); down = cbind(h_c, h_d); h_matrix = rbind(up, down); h_matrix 
  
  # Matrix variable
  h_y = rbind((t(X) %*% solve(R) %*% y), (t(Z) %*% solve(R) %*% y)); h_y
  
  # BLUE and BLUP
  all = solve(h_matrix) %*% h_y
  
  # Retornando os valores
  list(BLUE_BLUP = all)
}

#Função para calcular BLUP e BLUE - Henderson
bb_Hend2 <- function(y, X, Z, sigma2_e, sigma2_u) {
  #Dimensões das matrizes
  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(Z)
  
  # Matriz de variância do resíduo
  R <- diag(sigma2_e, n, n)
  
  # Matriz de variância de efeitos aleatórios
  G <- diag(sigma2_u, q, q)
  
  # Covariance matrix for the vector of observation y
  V_inv = solve(R) - (solve(R) %*% Z %*% solve((t(Z) %*% solve(R) %*% Z) + solve(G)) %*% t(Z)  %*% solve(R)) 
  
  # BLUE
  b = (solve(t(X) %*% V_inv %*% X)) %*% (t(X) %*% V_inv %*% y)
  
  # BLUP
  u = G %*% t(Z) %*% V_inv %*% (y - (X %*% b))
  
  # Retornando os valores
  list(BLUE = b, BLUP = u)
}

# Dados
y = matrix(c(9, 12, 11, 6, 7, 14),
           nrow = 6, ncol = 1, byrow = T); y

X = matrix(c(1, 0,
             0, 1,
             1, 0,
             1, 0,
             1, 0,
             0, 1),
           nrow = 6, ncol = 2, byrow = T); X

Z = matrix(c(1, 0, 0,
             1, 0, 0,
             0, 1, 0,
             0, 1, 0,
             0, 0, 1,
             0, 0, 1),
           nrow = 6, ncol = 3, byrow = T); Z

b = matrix()

u = matrix()

# Variâncias assumidas
sigma2_e <- 6  # Variância do resíduo
sigma2_u <- 8/4 # Variância dos efeitos aleatórios

#Chamando a função para calcular BLUP e BLUE
results.01 <- bb_MLE(y, X, Z, sigma2_e, sigma2_u); results.01
results.02 <- bb_Hend1(y, X, Z, sigma2_e, sigma2_u); results.02
results.03 <- bb_Hend2(y, X, Z, sigma2_e, sigma2_u); results.03

#Comprarando os tempos de execução
time.01 <- system.time(bb_MLE(y, X, Z, sigma2_e, sigma2_u)); time.01
time.02 <- system.time(bb_Hend1(y, X, Z, sigma2_e, sigma2_u)); time.02
time.03 <- system.time(bb_Hend2(y, X, Z, sigma2_e, sigma2_u)); time.03
