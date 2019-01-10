package org.intermine.bio.dataconversion;

/*
 * Copyright (C) 2002-2016 FlyMine
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  See the LICENSE file for more
 * information or http://www.gnu.org/copyleft/lesser.html.
 *
 */

import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.template.Sequence;
import org.intermine.model.bio.BioEntity;
import org.intermine.model.bio.DataSet;
import org.intermine.model.bio.Organism;
import org.intermine.objectstore.ObjectStoreException;

/**
 * A loader that works for FASTA files from TGI, used to extract
 * the FASTA header and store as briefDescription
 * @author Vivek Krishnakumar
 */
public class MTGIFastaLoaderTask extends FastaLoaderTask
{
   /**
     * {@inheritDoc}
     */
    @Override
    protected void extraProcessing(Sequence bioJavaSequence,
            org.intermine.model.bio.Sequence flymineSequence,
            BioEntity bioEntity, Organism organism, DataSet dataSet)
        throws ObjectStoreException {

            String header = ((DNASequence) bioJavaSequence).getOriginalHeader();

            // set the briefDescription attribute of a TC SequenceFeature
            // trim() the header to remove leading/trailing whitespaces
            bioEntity.setFieldValue("briefDescription", header.trim());
//        }
    }
}
