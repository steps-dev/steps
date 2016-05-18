# Check for OSGeo4W.bat binary

.onLoad <-
  function(libname, pkgname) {
    op <-
      base::options()
    op.dlmpr <-
      base::list(
        dlmpr.path = {
          if (.Platform$OS.type == "windows") {
            base::Sys.which("OSGeo4W.bat")
          } else {
            base::Sys.which("gdal_polygonize.py")
          }
        }
      )
    toset <-
      op.dlmpr   %>%
      base::names(.) %>%
      magrittr::is_in(
        x     = .,
        table = base::names(op)
      ) %>%
      magrittr::not(.);
    if (base::any(toset)) {
      base::options(op.dlmpr[toset]);
    };
    base::invisible();
  }
