
ADRminerServer <- function(what=c("gps")){
  what <- match.arg(what)
  if(what=="gps"){
    .gpsServer()
  }
  
  return(invisible())
}

###############
## .gpsServer
###############
## hidden function - DAPC server
.gpsServer <- function(){
  runApp(system.file("gpsServer",package="ADRminer"))
}
