### Read in nexus formatted output from splitstree.
read_splitstree_nexus = function(fname) {
    section = 'none'
    subsection = 'none'
    nexus_file_string = scan(fname, 'rb', sep='\n', what='character')
    i = 0
    line = ''
    while(i < length(nexus_file_string)) {
        i = i + 1
        line = trimws(nexus_file_string[i], 'both')
        if(length(line) == 0) {
            next
        } else if(grepl('^END;', line)) {
            section = 'none'
            next
        } else if(grepl('^BEGIN', line)) {
            section = gsub('BEGIN (.+);$', '\\1', line)
            next
        }
        if(section == 'Taxa') {
            if(grepl('^DIMENSIONS', line)) {
                taxa = character(as.numeric(gsub('.+ntax=(.+);', '\\1', line)))
            } else if(grepl('^TAXLABELS|;$', line)) {
                next
            } else {
                t_number = as.numeric(gsub('\\[([0-9]+)\\].+', '\\1', line))
                t_name = gsub("^'|'$", '', gsub('\\[[0-9]+\\] ', '', line))
                taxa[t_number] = t_name
            }
        } else if(section == 'Distances') {
            ## Don't do anything with this for now
            next
        } else if(section == 'Splits') {
            ## Don't do anything with this for now
            next
        } else if(section == 'Network') {
            if(grepl('^DIMENSIONS', line)) {
                n_taxa_plot = as.numeric(gsub('.+ntax=([0-9]+).+', '\\1', line))
                n_vertex = as.numeric(gsub('.+nvertices=([0-9]+).+', '\\1', line))
                n_edge = as.numeric(gsub('.+nedges=([0-9]+).+', '\\1', line))
                edges = data.frame('number'=numeric(n_edge),
                    'v1'=numeric(n_edge),
                    'v2'=numeric(n_edge),
                    'w'=numeric(n_edge),
                    stringsAsFactors=FALSE)
                vertices = data.frame(label=character(n_vertex),
                    number=numeric(n_vertex),
                    leaf=logical(n_vertex),
                    stringsAsFactors=FALSE)
                vertex_data = data.frame(number=numeric(n_vertex),
                    x=numeric(n_vertex),
                    y=numeric(n_vertex))
                subsection = 'none'
            } else if(line == ';') {
                subsection = 'none'
                next
            } else {
                if(line == 'TRANSLATE') {
                    subsection = 'translate'
                    next
                } else if(line == 'EDGES') {
                    subsection = 'edges'
                    next
                } else if(line == 'VERTICES') {
                    subsection = 'vertices'
                    next
                }
                if(subsection == 'translate') {
                    v_number = as.numeric(unlist(strsplit(line, ' '))[1])
                    vertices[v_number, 'label'] = gsub("',$", '', gsub("[0-9]+ '", '', line))
                } else if(subsection == 'edges') {
                    edge_data = unlist(strsplit(line, ' '))
                    e_number = as.numeric(edge_data[1])
                    v1 = as.numeric(edge_data[2])
                    v2 = as.numeric(edge_data[3])
                    w = as.numeric(gsub('w=', '',
                        edge_data[grepl('^w=', edge_data)]))
                    edges[e_number, ] = c(e_number, v1, v2, w)
                } else if(subsection == 'vertices') {
                    v_data = unlist(strsplit(line, ' '))
                    v_number = as.numeric(v_data[1])
                    vertex_data[v_number, 'x'] = as.numeric(v_data[2])
                    vertex_data[v_number, 'y'] = as.numeric(v_data[3])
                    vertices[v_number, 'number'] = v_number
                    vertex_data[v_number, 'number'] = v_number
                }
            }
        }
    }

    #parents = unique(edges[, 1])
    #children = unique(edges[, 2])
    #leaves = children[!children %in% parents]
    tab = table(c(edges$v1, edges$v2))
    leaves = names(tab[tab == 1])
    vertices$leaf[vertices$number %in% leaves] = TRUE
    list('vertices'=vertices, 'edges'=edges, 'vertex_data'=vertex_data)
}

add_vertex_data = function(net_obj, dat, match_label=FALSE) {
    numbers = net_obj[['vertices']][, 'number']
    labels = net_obj[['vertices']][, 'label']
    if(!match_label) {
        stopifnot('number' %in% names(dat))
        mcol = 'number'
        m = numbers
    } else {
        stopifnot('label' %in% names(dat))
        mcol = 'label'
        m = labels
    }

    old_data = net_obj[['vertex_data']]
    old_data = old_data[match(numbers, old_data[, 'number']), ]
    old_names = names(old_data)
    old_names = old_names[old_names != 'number']
    new_names = names(dat)
    new_names = new_names[!(new_names %in% c('number', 'label'))]
    all_names = unique(c(old_names, new_names))
    new_data = structure(vector('list', length=length(all_names)+1),
                         names=c('number', all_names))
    
    new_data[['number']] = numbers
    for(nm in old_names) {
        new_data[[nm]] = old_data[[nm]]
    }
    for(nm in new_names) {
        new_data[[nm]] = dat[match(m, dat[, mcol]), nm]
    }


    net_obj[['vertex_data']] = do.call(data.frame,
                                       c(new_data,
                                         stringsAsFactors=FALSE,
                                         check.names=FALSE))
    net_obj
}


## Works!
plot_network = function(net_obj,
    label=NULL,
    pts=NULL,
    pts_lvs_only=TRUE,
    symbol=NULL,
    symbol_lvs_only=TRUE,
    col.lab='black',
    cex.lab=0.3,
    font.lab=1,
    adj.lab=NULL,
    pos.lab=NULL,
    offset.lab=0.5,
    vfont.lab=NULL,
    col.edge=par('fg'),
    lwd.edge=1,
    lty.edge=1,
    ...) {

    vertices = net_obj$vertices
    vdata = net_obj$vertex_data
    stopifnot(all(vdata$number == vertices$number))
    edges = net_obj$edges
    x = vdata$x
    y = vdata$y

    plot(0, 0, xlim=c(min(x), max(x)),
         ylim=c(min(y), max(y)),
         type='n', bty='n', xaxt='n', yaxt='n', mar=c(0, 0, 0, 0),
         xlab='', ylab='', ...)

    for(i in 1:nrow(edges)) {
        v1 = edges[i, 'v1']
        v2 = edges[i, 'v2']
        x1 = x[vertices$number == v1]
        y1 = y[vertices$number == v1]
        x2 = x[vertices$number == v2]
        y2 = y[vertices$number == v2]
        segments(x1, y1, x2, y2,
                 lwd=lwd.edge, lty=lty.edge, col=col.edge)
    }

    if(!is.null(label)) {
        if(is.character(label)) {
            labs = label
        } else {
            labs = vertices[, 'label']
        }
        text(vdata[, 'x'], vdata[, 'y'],
             labs, cex=cex.lab, col=col.lab,
             font=font.lab, adj=adj.lab, pos=pos.lab,
             offset=offset.lab, vfont=vfont.lab)
    }

    if(!is.null(pts)) {
        xx = x
        yy = y
        if(pts_lvs_only) {
            xx = x[vertices[, 'leaf']]
            yy = y[vertices[, 'leaf']]
        }
        do.call(points, c(list('x'=x, 'y'=y), pts))
    }

    if(!is.null(symbol)) {
        xx = x
        yy = y
        if(symbol_lvs_only) {
            xx = x[vertices[, 'leaf']]
            yy = y[vertices[, 'leaf']]
        }
        do.call(symbols, c(list('x'=xx, 'y'=yy, add=TRUE), symbol))
    }

}
