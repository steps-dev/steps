# Check for OSGeo4W.bat and gdal_polygonize.py
.onLoad <-
  function(libname, pkgname) {
    op <-
      base::options()
    op.dlmpr <-
      base::list(
        cmd.path = {
          if (.Platform$OS.type == "windows") {
            base::Sys.which("OSGeo4W.bat")
          } else {
            'python'
          }
        },
        py.path = {
          system.file("python", "gdal_polygonize.py", package = "dlmpr");
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
