import React, { useState, useEffect } from 'react';
import { makeStyles } from '@material-ui/core/styles';
import Paper from '@material-ui/core/Paper';
import Tabs from '@material-ui/core/Tabs';
import Tab from '@material-ui/core/Tab';
import DropdownButton from 'react-bootstrap/DropdownButton';
import Dropdown from 'react-bootstrap/Dropdown';
import BlastReportDescription from './reportDescription.jsx';
import BlastReportAlignments from './reportAlignments.jsx';

const useStyles = makeStyles(theme => ({
    root: {
        width: '100%',
        marginTop: '20px',
        overflowX: 'auto',
        position: 'sticky',
        top: '0px',
    },
}));

function BlastReport({ report }) {
    const [tab, setTab] = useState(0);
    const [algnItem, setAlgnItem] = useState(null);

    function handleTabClick(event, newValue) {
        setTab(newValue);
    }

    useEffect(() => {
        if (algnItem !== null && tab === 1) {
            $('html, body').animate(
                {
                    scrollTop: $('#dln_' + algnItem).offset().top - 60,
                },
                300
            );
            setAlgnItem(null);
        }
    }, [tab]);

    function handleItemClick(event) {
        event.preventDefault();
        setTab(1);
        setAlgnItem(event.target.getAttribute('href'));
    }

    const classes = useStyles();

    return (
        <div>
            <Paper square className={classes.root}>
                <Tabs
                    value={tab}
                    onChange={handleTabClick}
                    indicatorColor="primary"
                    textColor="primary"
                >
                    <Tab label="Descriptions" />
                    <Tab label="Alignments" />
                </Tabs>
            </Paper>
            {tab === 0 && (
                <div>
                    <BlastReportDescription
                        report={report.report.results}
                        onClick={handleItemClick}
                    />
                </div>
            )}
            {tab === 1 && (
                <div>
                    <BlastReportAlignments report={report.report.results} />
                </div>
            )}
        </div>
    );
}

export default BlastReport;
