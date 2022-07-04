# GENERATE SIGNIFICANCE INDEX

  defast = function(x) {
    if (x > 0.1) {
      ast = " "
    } else
    {
      if (x > 0.05) {
        ast = "."
      } else
      {
        if (x > 0.01) {
          ast = "*"
        } else
        {
          if (x > 0.001) {
            ast = "**"
          } else
          {
            {
              ast = "***"
            }
          }
        }
      }
    }
    return(ast)
  }

