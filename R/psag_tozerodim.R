psag_tozerodim <- function(X, V, type = '2'){
  d = nrow(X) - 1
  n = ncol(X)

  m = rowMeans(X) # backwards mean
  idx = 2
  removed_vertex = V[,idx]
  m_y = mat_inverse(V) %*% m
  a = m_y[2]/m_y[1]

  Vhat = matrix(m, d+1, 1)
  rownames(Vhat) = rownames(X)
  colnames(Vhat) = 'V1'
  Xhat = Vhat %*% matrix(1, 1, n)
  Res = X - Xhat

  Y = mat_inverse(V) %*% X
  Xhat_reduced = matrix(1, nrow = 1, ncol = n, dimnames = list('V1', NULL))
  Yhat = m_y %*% matrix(1, 1, n)
  Res_y = Y - Yhat
  scores = sign(Res_y[nrow(Res_y),]) * sqrt(colSums(Res^2))
  RSS = sum(scores^2)

  return(list(Xhat = Xhat, Vhat = Vhat, Xhat_reduced = Xhat_reduced,
              idx = idx, a = a, removed_vertex = removed_vertex,
              residuals = Res, scores = scores, RSS = RSS))
}
