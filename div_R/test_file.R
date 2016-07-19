ge <- c("a", "a", "a", "a",
        "b", "b", "b", "b",
        "c", "c", "c", "c")

sirnas <- c("a1", "a1", "a2", "a2",
            "b1", "b1", "b2", "b2",
            "c1", "c1", "c2", "c2")
re <- rnorm(12)

level <- "gene"

.hypertest(ge, sirnas, re, "gene")
