if( !isGeneric("result") )
    setGeneric("result", function(object) standardGeneric("result"))

if( !isGeneric("p.value") )
    setGeneric("p.value", function(gt) standardGeneric("p.value"))

if( !isGeneric("fit") )
    setGeneric("fit", function(x) standardGeneric("fit"))

if( !isGeneric("names") ) setGeneric("names")

if( !isGeneric("names<-") ) setGeneric("names<-")

if( !isGeneric("sort") ) setGeneric("sort")

if( !isGeneric("z.score") )
    setGeneric("z.score", function(x) standardGeneric("z.score"))

if( !isGeneric("scale")) setGeneric("scale")

setGeneric("combine")

if( !isGeneric("plot") )
    setGeneric("plot", function(x, y, ...) standardGeneric("plot"))
