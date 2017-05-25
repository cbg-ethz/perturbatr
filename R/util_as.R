# knockout: analysis of high-throughput gene perturbation screens
#
# Copyright (C) 2015 - 2016 Simon Dirmeier
#
# This file is part of knockout
#
# knockout is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# knockout is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with knockout. If not, see <http://www.gnu.org/licenses/>.


setAs(
  "data.table",
  "knockout.data",
  function(from)
  {
    req.els <- list("GeneSymbol" = "character",
                    "Control"    = "integer",
                    "siRNAIDs"   = "character",
                    "Readout"    = c("numeric", "integer"))

    col.names <- colnames(from)
    names <- names(req.els)
    for (i in seq(req.els))
    {
      name <- names[i]
      if (!name %in% col.names)
      {
        stop(paste0("object misses a column called: ",  "'", name, "'."))
      }
      if (!class(from[[name]]) %in% req.els[[i]])
      {
        stop(paste0("column '", name, "'", " has not class ", req.els[[i]], "."))
      }
    }

    els <- list("Virus"        = NA_character_,
                "Replicate"    = NA_integer_,
                "Plate"        = NA_integer_,
                "RowIdx"       = NA_integer_,
                "ColIdx"       = NA_integer_,
                "GeneSymbol"   = NA_character_,
                "ReadoutType"  = NA_character_,
                "Control"      = NA_integer_,
                "Library"      = NA_character_,
                "siRNAIDs"     = NA_character_,
                "Screen"       = NA_character_,
                "Cell"         = NA_character_,
                "ScreenType"   = NA_character_,
                "Design"       = NA_character_,
                "Entrez"       = NA_character_,
                "Readout"      = NA_real_,
                "ReadoutClass" = NA_character_,
                "NumCells"     = NA_character_)
    names <- names(els)
    for (i in seq(els))
    {
      name <- names[i]
      if (!(name %in% col.names))
      {
        from[[name]] <- els[[i]]
      }
    }
    return(methods::new("knockout.raw.data", .data=from))
})

