"""File output for field(s) value on a grid.
"""


class XDMFWriter(object):

    def __init__(self, h5filename, dimension, resolution, origin, space_step,
                 dataset_names, ite, time):
        """
        Parameters
        ----------
        h5filename : nom du fichier h5 contenant les donnees
        dimension du domaine
        resolution de la grille
        origin : coordonnees de l'origine du domain
        space_step : pas d'espace
        dataset_names : liste des datasets qu'on souhaite ecrire dans le xdmf
        ite : numero de l'iteration courante
        time : temps

        """
        self.xmffilename = h5filename.split('.')[0] + '.xmf'
        res = list(resolution)

        f = open(self.xmffilename, 'w')
        f.write("<?xml version=\"1.0\" ?>\n")
        f.write("<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\">\n")
        f.write("<Xdmf Version=\"2.0\">\n")
        f.write(" <Domain>\n")
        f.write("  <Grid Name=\"CellTime\" GridType=\"Collection\" ")
        f.write("CollectionType=\"Temporal\">\n")
        f.write(self._write_grid_attributes(dimension, res,
                                            origin, space_step,
                                            dataset_names, ite,
                                            time, h5filename))
        f.write("  </Grid>\n")
        f.write(" </Domain>\n")
        f.write("</Xdmf>\n")
        print('Ecriture du fichier ' + self.xmffilename)
        f.close()

    def _write_grid_attributes(self, dimension, resolution, origin, space_step,
                               dataset_names, ite,
                               time, filename):
        """
        Write XDMF header into a file

        Returns:
        --------
        string
            the xml-like header.

        """
        assert isinstance(resolution, list)
        assert isinstance(origin, list)
        assert isinstance(space_step, list)
        # The header (xml-like), saved in a string.
        xml_grid = ""
        if dimension == 2:
            topo_type = "2DCORECTMesh"
            geo_type = "ORIGIN_DXDY"
        elif dimension == 3:
            topo_type = "3DCORECTMesh"
            geo_type = "ORIGIN_DXDYDZ"
        xml_grid += "   <Grid Name=\"Iteration {0:03d}\"".format(ite)
        xml_grid += " GridType=\"Uniform\">\n"
        xml_grid += "    <Time Value=\"{0}\" />\n".format(time)
        xml_grid += "    <Topology TopologyType=\"" + str(topo_type) + "\""
        xml_grid += " NumberOfElements=\""
        resolution.reverse()
        origin.reverse()
        xml_grid += XDMFWriter._list_format(resolution) + " \"/>\n"
        xml_grid += "    <Geometry GeometryType=\"" + geo_type + "\">\n"
        xml_grid += "     <DataItem Dimensions=\"" + str(dimension) + " \""
        xml_grid += " NumberType=\"Float\" Precision=\"8\" Format=\"XML\">\n"
        xml_grid += "     " + XDMFWriter._list_format(origin) + "\n"
        xml_grid += "     </DataItem>\n"
        xml_grid += "     <DataItem Dimensions=\"" + str(dimension) + " \""
        xml_grid += " NumberType=\"Float\" Precision=\"8\" Format=\"XML\">\n"
        step = space_step
        step.reverse()
        xml_grid += "     " + XDMFWriter._list_format(step) + "\n"
        xml_grid += "     </DataItem>\n"
        xml_grid += "    </Geometry>\n"
        # Append dataset parameters
        for name in dataset_names:
            xml_grid += "    <Attribute Name=\""
            xml_grid += name + "\""
            xml_grid += " AttributeType=\"Scalar\" Center=\"Node\">\n"
            xml_grid += "     <DataItem Dimensions=\""
            xml_grid += XDMFWriter._list_format(resolution) + " \""
            xml_grid += " NumberType=\"Float\" Precision=\"8\" Format=\"HDF\""
            xml_grid += " Compression=\"Raw\">\n"  #
            xml_grid += "      " + filename.split('/')[-1]
            xml_grid += ":/" + name
            xml_grid += "\n     </DataItem>\n"
            xml_grid += "    </Attribute>\n"
        xml_grid += "   </Grid>\n"
        return xml_grid

    @staticmethod
    def _list_format(l):
        """Format a list to the xml output.
        Removes the '[]()' and replace ',' with ' ' in default str.

        Parameters
        ----------
        l : list to format


        """
        buff = str(l).replace(',', ' ').replace('[', '')
        return buff.replace(']', '').replace('(', '').replace(')', '')
