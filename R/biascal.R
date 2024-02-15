biascal_new <- function(invcor, Amat, dir, val, k, p, Q){
  aprime <- invcor

  ## Max
  mycop_max <- cop(f=linfun(a=aprime, name="Obj.function"),
                   lc=lincon(A=Amat, dir=dir, val=val, name=seq(1,(k+p)*2)),
                   qc=quadcon(Q=Q, d=0, dir="<=", val=1, use=TRUE),
                   max=TRUE)
  ## Min
  mycop_min <- cop(f=linfun(a=aprime, name="Obj.function"),
                   lc=lincon(A=Amat, dir=dir, val=val, name=seq(1,(k+p)*2)),
                   qc=quadcon(Q=Q, d=0, dir="<=", val=1, use=TRUE),
                   max=FALSE)

  res_max <- solvecop(op=mycop_max, solver="alabama", quiet=TRUE)
  res_min <- solvecop(op=mycop_min, solver="alabama", quiet=TRUE)
  maxbias <- validate(op=mycop_max, sol=res_max, quiet=TRUE)$obj.fun
  minbias <- validate(op=mycop_min, sol=res_min, quiet=TRUE)$obj.fun

  return(c(minbias, maxbias))
}

biascal_new_forpdp <- function(data, x, pdvalue_i, i, k, p, expindx, bound){
  means <- colMeans(x)
  means[(k+expindx)] <- pdvalue_i
  invcor <- as.matrix(solve(cor(data[,-1])))
  linconst <- invcor %*% means

  val <- bound
  dir <- c(rep("<=",(k+p)), rep(">=",(k+p)))
  Amat <- as.matrix(rbind(diag(1,k+p), diag(1,k+p)))

  myQ3 <- invcor

  mycop_max <- cop(f=linfun(a=linconst, name="Obj.function"),
                   lc=lincon(A=Amat, dir=dir, val=val, name=seq(1,(k+p)*2)),
                   qc=quadcon(Q=myQ3, d=0, dir="<=", val=1, use=TRUE),
                   max=TRUE)

  mycop_min <- cop(f=linfun(a=linconst, name="Obj.function"),
                   lc=lincon(A=Amat, dir=dir, val=val, name=seq(1,(k+p)*2)),
                   qc=quadcon(Q=myQ3, d=0, dir="<=", val=1, use=TRUE),
                   max=FALSE)

  res_max <- solvecop(op=mycop_max, solver="alabama", quiet=TRUE)
  res_min <- solvecop(op=mycop_min, solver="alabama", quiet=TRUE)
  maxbias <- validate(op=mycop_max, sol=res_max, quiet=TRUE)$obj.fun
  minbias <- validate(op=mycop_min, sol=res_min, quiet=TRUE)$obj.fun

  return(c(minbias, maxbias))
}
