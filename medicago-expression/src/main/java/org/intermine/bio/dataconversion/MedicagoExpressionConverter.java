package org.intermine.bio.dataconversion;

/*
 * Copyright (C) 2002-2011 FlyMine
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  See the LICENSE file for more
 * information or http://www.gnu.org/copyleft/lesser.html.
 *
 */
import java.io.File;
import java.io.FileReader;

import java.io.Reader;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.math.BigDecimal;

import org.intermine.dataconversion.ItemWriter;
import org.intermine.metadata.Model;
import org.intermine.xml.full.Item;
import org.intermine.objectstore.ObjectStoreException;
import org.intermine.util.FormattedTextParser;
import org.apache.commons.lang.StringUtils;
import org.apache.log4j.Logger;
import org.apache.tools.ant.BuildException;


/**
 *
 * @author
 */
public class MedicagoExpressionConverter extends BioFileConverter
{
    //
    private static final String DATASET_TITLE = "RNA-seq expression";
    private static final String DATA_SOURCE_NAME = "MTGD";
    private static final Logger LOG = Logger.getLogger(MedicagoExpressionConverter.class);
    private Item organism;
    private static final String TAXON_ID = "3880";
    private Map<String, String> mrnas = new HashMap<String, String>();
    private Map<String, String> terms = new HashMap<String, String>();


    /**
     * Constructor
     * @param writer the ItemWriter used to handle the resultant items
     * @param model the Model
     */
    public MedicagoExpressionConverter(ItemWriter writer, Model model) {
        super(writer, model, DATA_SOURCE_NAME, DATASET_TITLE);
        organism = createItem("Organism");
        organism.setAttribute("taxonId", TAXON_ID);
        try {
            store(organism);
        } catch (ObjectStoreException e) {
            throw new RuntimeException(e);
        }
    }
    public void process (Reader reader)  throws ObjectStoreException {
        Iterator<?> tsvIter;
        try {
            tsvIter = FormattedTextParser.parseTabDelimitedReader(reader);
        } catch (Exception e) {
            throw new BuildException("cannot parse file: " + getCurrentFile(), e);
        }
        while (tsvIter.hasNext()) {
            String[] line = (String[]) tsvIter.next();
            if(line.length < 3){
                continue;
            }
            String id = line[0];	// Medtr1g005000
            Item result = createItem("RNASeqResult");
	    BigDecimal bd = new BigDecimal(line[1]); 	// 0.130758
	    bd = bd.setScale(2,BigDecimal.ROUND_HALF_UP);
            String score = bd.toString() ;
            String stage = line[2];
            if (StringUtils.isNotEmpty(score)) {
		 try {
		     //    Float score = Float.valueOf(fpkm).floatValue();
                    result.setAttribute("expressionLevel", score);
		     } catch (NumberFormatException e) {
		     LOG.warn("bad score: " + score, e);

		     }
            }
             if (StringUtils.isNotEmpty(stage)) {
                try {
                    result.setAttribute("stage", stage);
                } catch (NumberFormatException e) {
                    LOG.warn("bad stage: " + stage, e);
                }
            }
            String mrna = getMRNA(id);
            if (StringUtils.isNotEmpty(mrna)) {
                //result.setReference("transcript", mrna);
		result.setReference("mrna", mrna);
                store(result);
            }
        }
    }
    private String getMRNA(String fbgn) throws ObjectStoreException {
        if (StringUtils.isEmpty(fbgn)) {
            return null;
        }

	if(mrnas.containsKey(fbgn)) {
	    return mrnas.get(fbgn);
	}
	//Item mrna = createItem("Transcript");
        Item mrna = createItem("MRNA");
	mrna.setAttribute("primaryIdentifier", fbgn);
	mrna.setReference("organism", organism);

	String refId = mrna.getIdentifier();
	mrnas.put(fbgn, refId);
	store(mrna);
        return refId;
    }

    class Stage {
        protected String name;
        protected String category;

        public Stage(String name, String category) {
            this.name = name;
            this.category = category;
        }
    }
}
