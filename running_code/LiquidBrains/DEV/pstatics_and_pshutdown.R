load('~/research/gits/AnTracks/results/pstatic.RData')
load('~/research/gits/AnTracks/results/pshutdown.RData')

ps <- rbindlist(list(shutdown = pshutdown, static = pstatic), idcol = TRUE)
mltps <- melt(ps, id.vars = c('.id', 'rho'))

grid.arrange(
ggplot(data = mltps, aes(rho, value, color = factor(.id))) + 
	geom_point() + #geom_line()+
	facet_wrap(~ factor(variable))+
	geom_smooth(formula = y ~ x, method = 'lm', se = FALSE) +
	scale_color_manual(values = c('mediumpurple', 'gold3')),

ggplot(data = ps, aes(rho, pnest / pfood, color = factor(.id))) + 
	geom_point() + #geom_line()+
	geom_smooth(formula = y ~ x, method = 'lm', se = FALSE) +
	scale_color_manual(values = c('mediumpurple', 'gold3')),
ggplot(data = ps, aes(rho, not_found / pnest, color = factor(.id))) + 
	geom_point() + #geom_line()+
	geom_smooth(formula = y ~ x, method = 'lm', se = FALSE) +
	scale_color_manual(values = c('mediumpurple', 'gold3')),

ggplot(data = ps, aes(rho, not_found / pfood, color = factor(.id))) + 
	geom_point() + #geom_line()+
	geom_smooth(formula = y ~ x, method = 'lm', se = FALSE) +
	scale_color_manual(values = c('mediumpurple', 'gold3'))
)

cor.test(ps[.id == 'static', rho], ps[.id == 'static', not_found])
cor.test(ps[.id == 'static', rho], ps[.id == 'static', not_found / pnest])
cor.test(ps[.id == 'static', rho], ps[.id == 'static', not_found / pfood])
cor.test(ps[.id == 'shutdown', rho], ps[.id == 'shutdown', pnest / pfood]) # R = -0.33, p = 0.14

cor.test(ps[.id == 'shutdown', rho], ps[.id == 'shutdown', pnest])
cor.test(ps[.id == 'shutdown', rho], ps[.id == 'shutdown', pfood]) # R = 0.35, p = 0.12
cor.test(ps[.id == 'static', rho], ps[.id == 'static', pnest])
cor.test(ps[.id == 'static', rho], ps[.id == 'static', pfood])


#################################################################
##############################     ##############################
##############################     ##############################
#################################################################


load('~/research/gits/AnTracks/results/pstatic.RData')
load('~/research/gits/AnTracks/results/pshutdown.RData')

ps <- rbindlist(list(shutdown = pshutdown, static = pstatic), idcol = TRUE)
mltps <- melt(ps, id.vars = c('.id', 'rho'))

grid.arrange(
	ggplot(data = mltps, aes(rho, value, color = factor(.id))) + 
		geom_point() + #geom_line()+
		facet_wrap(~ factor(variable))+
		geom_smooth(formula = y ~ x, method = 'lm', se = FALSE) +
		scale_color_manual(values = c('mediumpurple', 'gold3')),
	
	ggplot(data = ps, aes(rho, pnest / pfood, color = factor(.id))) + 
		geom_point() + #geom_line()+
		geom_smooth(formula = y ~ x, method = 'lm', se = FALSE) +
		scale_color_manual(values = c('mediumpurple', 'gold3')),
	ggplot(data = ps, aes(rho, not_found / pnest, color = factor(.id))) + 
		geom_point() + #geom_line()+
		geom_smooth(formula = y ~ x, method = 'lm', se = FALSE) +
		scale_color_manual(values = c('mediumpurple', 'gold3')),
	
	ggplot(data = ps, aes(rho, not_found / pfood, color = factor(.id))) + 
		geom_point() + #geom_line()+
		geom_smooth(formula = y ~ x, method = 'lm', se = FALSE) +
		scale_color_manual(values = c('mediumpurple', 'gold3'))
)

cor.test(ps[.id == 'static', rho], ps[.id == 'static', not_found])
cor.test(ps[.id == 'static', rho], ps[.id == 'static', not_found / pnest])
cor.test(ps[.id == 'static', rho], ps[.id == 'static', not_found / pfood])
cor.test(ps[.id == 'shutdown', rho], ps[.id == 'shutdown', pnest / pfood]) # R = -0.33, p = 0.14

cor.test(ps[.id == 'shutdown', rho], ps[.id == 'shutdown', pnest])
cor.test(ps[.id == 'shutdown', rho], ps[.id == 'shutdown', pfood]) # R = 0.35, p = 0.12
cor.test(ps[.id == 'static', rho], ps[.id == 'static', pnest])
cor.test(ps[.id == 'static', rho], ps[.id == 'static', pfood])
