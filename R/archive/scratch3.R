plotmat2 <- 
  nih_plotmat |>
  pivot_longer(! rank) |>
  mutate(
    color = case_when(
      name == "pl_pl" ~ "#1f78b4",
      name == "npl_pl" ~ "#33a02c",
      name == "npl_npl" ~ "#b2df8a",
      name == "pl_npl" ~ "#a6cee3",
      name == "nih" ~ "#d95f02"
    ),
    lty = case_when(
      name == "pl_pl" ~ 2,
      name == "npl_pl" ~ 4,
      name == "npl_npl" ~ 2,
      name == "pl_npl" ~ 4,
      name == "nih" ~ 1
    )
  )

a1 <- plotmat2 |> 
  filter(name %in% c("npl_npl", "pl_npl")) |> 
  ggplot() + 
  geom_line(aes(x = rank, y = value, color = name), color = plotmat2$color, lty = plotmat2$lty, lwd = 1.25) + 
  ylim(1, 10000) + 
  scale_x_continuous(trans='log10') +
  scale_y_continuous(trans='log10') +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  ylab("") +
  xlab("Rank") + 
  ggtitle("Power Law Eta")
