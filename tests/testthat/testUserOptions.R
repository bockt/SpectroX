

testParseRange = function(){

  cat("--- testParseRange:  --- \n")

  stopifnot(parseRange("2:3") == c(2,3))

  stopifnot(parseRange("5") == c(5))

  stopifnot(parseRange("4:19") == c(4,19))

  cat("--- testParseRange:  PASS ALL TEST --- \n")

}



# run tests
testParseRange()
